const assert = require('node:assert/strict');
const test = require('node:test');
const fixtures = require('./fixtures/satellite-main.json');
const { loadPage } = require('./satellite-harness.cjs');

function close(actual, expected, tolerance = 1e-10) {
  assert.ok(Number.isFinite(actual), `Expected finite value, got ${actual}`);
  assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected}`);
}

for (const fixture of fixtures.cases) {
  test(`main regression: ${fixture.name} (${fixture.params.samples} samples)`, () => {
    const page = loadPage();
    const data = page.simulate(fixture.params);
    for (const [key, expected] of Object.entries(fixture.stats)) close(data[key], expected);
    for (const point of fixture.points) {
      for (const [key, expected] of Object.entries(point)) {
        if (key === 'index') continue;
        const actual = data[key][point.index];
        if (expected === null) assert.ok(Number.isNaN(actual));
        else if (typeof expected === 'boolean') assert.equal(actual, expected);
        else if (Array.isArray(expected)) expected.forEach((v, i) => close(actual[i], v));
        else close(actual, expected);
      }
    }
    let visible = 0;
    for (let i = 0; i < data.nframes; i++) {
      for (const key of ['times', 'subLat', 'subLon', 'azs', 'els']) {
        assert.ok(Number.isFinite(data[key][i]), `${key}[${i}]`);
      }
      assert.ok(data.subLat[i] >= -90 && data.subLat[i] <= 90);
      assert.ok(data.subLon[i] >= -180 && data.subLon[i] < 180);
      assert.ok(data.els[i] >= -90 && data.els[i] <= 90);
      assert.equal(data.visible[i], data.els[i] > 0);
      if (data.visible[i]) {
        visible++;
        close(Math.hypot(data.skyX[i], data.skyY[i]), 90 - data.els[i]);
      } else {
        assert.ok(Number.isNaN(data.skyX[i]) && Number.isNaN(data.skyY[i]));
      }
    }
    assert.equal(data.visibleFraction, visible / data.nframes);
  });
}

test('circular equatorial GEO is stationary at its nominal longitude', () => {
  const page = loadPage();
  for (const samples of [120, 720, 10000]) {
    const data = page.simulate({ ...page.params(), inclinationDeg: 0, samples });
    for (let i = 0; i < samples; i++) {
      close(data.subLat[i], 0);
      close(data.subLon[i], data.satLonDeg);
      close(data.els[i], data.els[0]);
    }
    assert.equal(data.visibleFraction, 1);
  }
});

test('inclined circular orbit agrees with analytic latitude and independent elevation geometry', () => {
  const page = loadPage();
  const params = { ...page.params(), samples: 10000, launchLatDeg: 75 };
  const data = page.simulate(params);
  const rad = Math.PI / 180;
  const phi = params.launchLatDeg * rad, lambda = params.launchLonDeg * rad;
  const earthE2 = (1 / 298.257223563) * (2 - 1 / 298.257223563);
  const N = 6378.137 / Math.sqrt(1 - earthE2 * Math.sin(phi) ** 2);
  const ox = N * Math.cos(phi) * Math.cos(lambda);
  const oy = N * Math.cos(phi) * Math.sin(lambda);
  const oz = N * (1 - earthE2) * Math.sin(phi);
  let visible = 0;
  for (let i = 0; i < data.nframes; i++) {
    const phase = 2 * Math.PI * i / (data.nframes - 1) + params.argPerigeeDeg * rad;
    close(data.subLat[i], Math.asin(Math.sin(params.inclinationDeg * rad) * Math.sin(phase)) / rad);
    const [x, y, z] = data.satPositions[i];
    const dx = x - ox, dy = y - oy, dz = z - oz;
    const up = dx * Math.cos(phi) * Math.cos(lambda) + dy * Math.cos(phi) * Math.sin(lambda) + dz * Math.sin(phi);
    const el = Math.asin(up / Math.hypot(dx, dy, dz)) / rad;
    close(data.els[i], el);
    if (el > 0) visible++;
  }
  assert.ok(data.visibleFraction > 0 && data.visibleFraction < 1);
  assert.equal(data.visibleFraction, visible / data.nframes);
});

test('rejects invalid sample counts before evaluating the orbit', () => {
  const page = loadPage();
  page.run('satEcef = () => { throw new Error("orbit should not be evaluated"); }');
  for (const samples of [NaN, Infinity, -Infinity, undefined, null, '', '720', 0, -1, 119, 720.5, 10000.1, 10001, 150000, 1e9]) {
    assert.throws(() => page.simulate({ ...page.params(), samples }), /Time resolution/);
  }
});

test('DOM input parsing, errors, recovery and restart preserve the previous valid simulation', () => {
  const page = loadPage();
  const samples = page.elements.get('samples');
  const error = page.elements.get('simulationError');
  const click = page.elements.get('restartBtn').listeners.click;
  const original = page.run('simData');
  const originalStats = page.elements.get('metricLat').textContent;
  for (const value of ['', ' ', 'NaN', 'Infinity', '1e309', '720.5', '720abc', '119', '10001', '150000']) {
    samples.value = value;
    click();
    assert.equal(page.run('simData'), original);
    assert.equal(page.elements.get('metricLat').textContent, originalStats);
    assert.equal(error.hidden, false);
    assert.match(error.textContent, /Time resolution/);
    assert.match(error.textContent, /previous simulation is unchanged/);
    assert.equal(samples.attributes['aria-invalid'], 'true');
    assert.equal(samples.focused, true);
  }
  for (const value of ['120', '121', '7.2e2', '10000']) {
    samples.value = value;
    click();
    assert.equal(page.run('simData.nframes'), Number(value));
    assert.equal(error.hidden, true);
    assert.equal(error.textContent, '');
    assert.equal(samples.attributes['aria-invalid'], undefined);
    assert.equal(page.run('frameIndex'), 0);
  }
  assert.equal(page.calls.alerts.length, 0);
  assert.match(page.html, /id="samples"[^>]*step="1"[^>]*min="120"[^>]*max="10000"[^>]*required/);
});

test('plot sampling preserves endpoints, extrema, visibility gaps and exact moving tips', () => {
  const page = loadPage();
  const data = page.simulate(fixtures.cases[2].params);
  const indices = Array.from(data.plot.indices);
  assert.equal(indices[0], 0);
  assert.equal(indices.at(-1), data.nframes - 1);
  assert.ok(indices.length < 1460, `unexpected plot size: ${indices.length}`);
  assert.equal(new Set(indices).size, indices.length);
  for (const field of ['subLat', 'subLon', 'els']) {
    let lo = 0, hi = 0;
    for (let i = 1; i < data.nframes; i++) {
      if (data[field][i] < data[field][lo]) lo = i;
      if (data[field][i] > data[field][hi]) hi = i;
    }
    assert.ok(indices.includes(lo) && indices.includes(hi));
  }
  let crossings = 0;
  for (let i = 1; i < data.nframes; i++) {
    if (data.visible[i] !== data.visible[i - 1]) {
      crossings++;
      assert.ok(indices.includes(i - 1) && indices.includes(i));
    }
  }
  assert.ok(crossings >= 2);
  const k = 3;
  assert.ok(!indices.includes(k));
  const tip = [123, 456];
  const trail = page.context.trailPoints(data.plot.ground, data.plot.indices, k, tip);
  assert.equal(trail.at(-1), tip);
  assert.equal(trail.length, indices.filter(i => i <= k).length + 1);
  for (const frame of [0, k, Math.floor(data.nframes / 2), data.nframes - 1]) {
    page.calls.lineTo = 0;
    page.draw(frame, data);
    assert.ok(page.calls.lineTo < 16000, `too much canvas work: ${page.calls.lineTo}`);
    assert.doesNotMatch(page.elements.get('metricTime').textContent, /NaN/);
  }
});

test('720-point plots retain every point and the animation still advances', () => {
  const page = loadPage();
  const data = page.run('simData');
  assert.equal(data.plot.indices.length, 720);
  data.plot.indices.forEach((value, i) => assert.equal(value, i));
  page.run('animate(100); animate(140); animate(180)');
  assert.equal(page.run('frameIndex'), 2);
});

test('high-eccentricity plotting retains longitude-wrap boundaries', () => {
  const page = loadPage();
  const data = page.simulate({ ...page.params(), samples: 10000,
    inclinationDeg: 60, eccentricity: 0.8, launchLatDeg: 65, argPerigeeDeg: 35 });
  const indices = new Set(data.plot.indices);
  let wraps = 0;
  for (let i = 1; i < data.nframes; i++) {
    if (Math.abs(data.subLon[i] - data.subLon[i - 1]) > 180) {
      wraps++;
      assert.ok(indices.has(i - 1) && indices.has(i));
    }
  }
  assert.ok(wraps > 0);
  // Compare the interpolated display line with every full-resolution sample.
  // This is geometric error in canvas pixels, not a rasterization test.
  let maxGroundError = 0;
  for (let j = 1; j < data.plot.indices.length; j++) {
    const a = data.plot.indices[j - 1], b = data.plot.indices[j];
    for (let i = a + 1; i < b; i++) {
      const f = (i - a) / (b - a);
      const lon = data.subLon[a] * (1 - f) + data.subLon[b] * f;
      const lat = data.subLat[a] * (1 - f) + data.subLat[b] * f;
      const error = Math.hypot((data.subLon[i] - lon) * 560 / (2 * data.lonMargin),
        (data.subLat[i] - lat) * 290 / (2 * data.latMargin));
      maxGroundError = Math.max(maxGroundError, error);
    }
  }
  assert.ok(maxGroundError < 0.5, `ground-track deviation: ${maxGroundError}px`);
});
