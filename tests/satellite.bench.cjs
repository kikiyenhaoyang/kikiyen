// Usage: node --expose-gc tests/satellite.bench.cjs [path/to/baseline/satellite.html]
const assert = require('node:assert/strict');
const v8 = require('node:v8');
const path = require('node:path');
const { performance } = require('node:perf_hooks');
const { loadPage } = require('./satellite-harness.cjs');

const median = values => values.sort((a, b) => a - b)[Math.floor(values.length / 2)];
const current = path.join(__dirname, '../Web/satellite.html');
const baseline = process.argv[2];
const rows = [];

for (const [version, filename] of [['before', baseline], ['after', current]]) {
  if (!filename) continue;
  const page = loadPage(filename);
  for (const samples of [720, 10000]) {
    const params = { ...page.params(), samples };
    for (let i = 0; i < 3; i++) page.simulate(params);
    const compute = [];
    for (let i = 0; i < 9; i++) {
      const start = performance.now();
      page.simulate(params);
      compute.push(performance.now() - start);
    }
    const data = page.simulate(params);
    const frame = Math.floor(samples / 2);
    for (let i = 0; i < 5; i++) page.draw(frame, data);
    const draw = [];
    for (let i = 0; i < 21; i++) {
      const start = performance.now();
      page.draw(frame, data);
      draw.push(performance.now() - start);
    }
    page.calls.lineTo = 0;
    page.draw(frame, data);
    let retainedKiB = null;
    if (global.gc) {
      // Hold multiple datasets to reduce noise from VM/JIT bookkeeping.
      const held = [];
      global.gc();
      const startHeap = v8.getHeapStatistics().used_heap_size;
      for (let i = 0; i < 10; i++) held.push(page.simulate(params));
      global.gc();
      retainedKiB = (v8.getHeapStatistics().used_heap_size - startHeap) / held.length / 1024;
      assert.equal(held[0].nframes, samples);
    }
    rows.push({ version, samples, computeMedianMs: median(compute),
      drawMedianMs: median(draw), lineToCalls: page.calls.lineTo, retainedKiB,
      plotPoints: data.plot?.indices.length ?? samples });
  }
}

if (baseline) {
  const before = loadPage(baseline), after = loadPage(current);
  for (const samples of [720, 10000]) {
    for (const overrides of [{}, { launchLatDeg: 65, inclinationDeg: 45,
      eccentricity: 0.3, argPerigeeDeg: 35, satLonDeg: 170, launchLonDeg: 100 }]) {
      const params = { ...before.params(), ...overrides, samples };
      const oldData = before.simulate(params), newData = after.simulate(params);
      // Compare every sample and statistic, excluding the new plotting cache.
      // NaN sky coordinates are intentional for below-horizon samples.
      for (const key of Object.keys(oldData)) {
        assert.equal(JSON.stringify(newData[key]), JSON.stringify(oldData[key]), `${samples}: ${key}`);
      }
    }
  }
  console.log('All simulation arrays and statistics match the baseline exactly (720/10,000; circular/eccentric).');
  try {
    before.simulate({ ...before.params(), samples: 150000 });
    console.log('This engine did not reproduce the baseline 150,000-point argument-limit error.');
  } catch (error) {
    console.log(`Baseline 150,000 samples: ${error.name}: ${error.message}`);
  }
  assert.throws(() => after.simulate({ ...after.params(), samples: 150000 }), /Time resolution/);
}
console.log(`Node ${process.version}; V8 ${process.versions.v8}; canvas calls are stubbed, not rasterized.`);
console.log(JSON.stringify(rows, null, 2));
