const test = require('node:test');
const assert = require('node:assert/strict');
const {loadSimulator, errors} = require('./grin-harness.cjs');

const base = {casimir: .030625, energy: 1.4517, alpha: .08, n0: 1.5, phase: 0};
const critical = {casimir: .98, energy: 1.30, alpha: .08, n0: 1.5, phase: 0};
function near(actual, expected, tolerance = 1e-12) {
  assert.ok(Math.abs(actual-expected) <= tolerance * Math.max(1,Math.abs(expected)),
    `${actual} differs from ${expected}`);
}
function physical(model, points = model.path) {
  assert.ok(points.length);
  for (const [X,Y,Z] of points) {
    assert.ok([X,Y,Z].every(Number.isFinite));
    assert.ok(X >= 0 && Y >= 0);
    near(X*Y-Z*Z, model.C);
    near(model.n0**2-model.alpha*X-Y, model.H**2);
    near(-model.alpha*X*X+model.b*X-model.C, Z*Z);
    assert.ok(X <= model.Xcut + 1e-12);
  }
  for (const X of model.turningPoints) {
    near(-model.alpha*X*X+model.b*X-model.C, 0);
    assert.ok(X >= 0);
    assert.ok(model.b-model.alpha*X >= -1e-12);
  }
}
function allGeometry(sim, model) {
  physical(model);
  physical(model, sim.run('energyCurvePoints(model,1)'));
  physical(model, sim.run('energyCurvePoints(model,-1)'));
  for (const sample of model.radiusSeries) {
    physical(model, [[sample.X,sample.Y,sample.Z]]);
    near(sample.r**2, sample.X);
  }
}

test('Guided: invariants, roots, phase, closed orbit and physical z period', () => {
  const sim = loadSimulator();
  for (const phase of [0, 37, 180, 359]) {
    const m = sim.configure({...base, phase});
    assert.equal(m.regime, 'guided');
    assert.equal(m.path.length, 420);
    assert.equal(m.turningPoints.length, 2);
    assert.ok(m.canAnimate);
    allGeometry(sim, m);
    near(m.zExtent, Math.PI*m.H/Math.sqrt(m.alpha));
    m.path[0].forEach((x,i) => near(x,m.path.at(-1)[i]));
    const curve = sim.run('energyCurvePoints(model)');
    near(curve[0][0],m.Xhi);
    near(curve.at(-1)[0],m.Xlo);
    near(curve[0][2],0);
    near(curve.at(-1)[2],0);
  }
  // Independent finite differences check the sign and z scale of all three ODEs.
  const m=sim.model(), i=120, h=1e-5, z=m.radiusSeries[i].z;
  const at = value => sim.run(`stateAtAngle(model,2*Math.sqrt(model.alpha)/model.H*${value}+model.phase0)`);
  const before=at(z-h), after=at(z+h), [X,Y,Z]=m.path[i];
  const derivative=after.map((x,j)=>(x-before[j])/(2*h));
  near(derivative[0],-2*Z/m.H,1e-8);
  near(derivative[1],2*m.alpha*Z/m.H,1e-8);
  near(derivative[2],(m.alpha*X-Y)/m.H,1e-8);
});

test('Critical report: one valid state, one double root, constant radius', () => {
  const sim=loadSimulator(), m=sim.configure(critical);
  assert.equal(m.regime,'critical');
  assert.equal(m.path.length,1);
  assert.equal(m.turningPoints.length,1);
  assert.equal(m.canAnimate,false);
  near(m.path[0][0],3.5);
  near(m.path[0][1],.28);
  assert.equal(m.path[0][2],0);
  near(m.b-2*m.alpha*m.Xlo,0);
  near(m.radiusSeries[0].r,m.radiusSeries.at(-1).r);
  allGeometry(sim,m);
});

test('Critical preset remains a double root after arbitrary knob settings', () => {
  const sim=loadSimulator();
  for (const n0 of [1.05,1.5,1.8]) for (const alpha of [.01,.08,.22]) {
    sim.configure({n0,alpha});
    sim.element('preset-critical').dispatch('click');
    const m=sim.model();
    assert.equal(m.regime,'critical');
    assert.equal(m.path.length,1);
    assert.equal(m.n0,1.5);
    assert.equal(m.alpha,.08);
    near(m.C,m.b*m.b/(4*m.alpha));
    assert.ok(sim.element('togglePlay').disabled);
    allGeometry(sim,m);
  }
});

test('Near-critical signs do not create a broad false critical band', () => {
  const sim=loadSimulator(), b=1.5**2-1.3**2, c=b*b/(4*.08);
  for (const offset of [1e-12,1e-8,1e-5]) {
    const guided=sim.configure({...critical,casimir:c-offset});
    assert.equal(guided.regime,'guided');
    allGeometry(sim,guided);
    const empty=sim.configure({...critical,casimir:c+offset});
    assert.equal(empty.regime,'no real trajectory');
    assert.equal(empty.path.length,0);
  }
  // A genuinely tiny positive C is not silently replaced by the C=0 cone.
  const small=sim.configure({...base,casimir:1e-10});
  assert.equal(small.regime,'guided');
  assert.ok(small.Xlo>0);
  allGeometry(sim,small);
});

test('Singular cone: both turning points, both Z signs and the energy constraint', () => {
  const sim=loadSimulator(), m=sim.configure({...base,casimir:0,energy:1.42});
  assert.equal(m.regime,'singular cone');
  assert.equal(m.path.length,420);
  assert.equal(m.Xlo,0);
  near(m.Xhi,m.b/m.alpha);
  assert.ok(m.path.some(p=>p[2]>0) && m.path.some(p=>p[2]<0));
  allGeometry(sim,m);
  const axis=sim.run('energyCurvePoints(model).at(-1)');
  near(axis[0],0); near(axis[1],m.b); near(axis[2],0);
  assert.ok(axis[1]>0); // The axis crossing is not the singular vertex.
});

test('Singular vertex only at C=0, H=n0; H>n0 has no physical state', () => {
  const sim=loadSimulator(), m=sim.configure({...base,casimir:0,energy:1.5});
  assert.equal(m.path.length,1);
  assert.equal(m.canAnimate,false);
  assert.deepEqual(Array.from(m.path[0]),[0,0,0]);
  allGeometry(sim,m);
  for (const casimir of [0,.8]) {
    const empty=sim.configure({...base,casimir,energy:1.58});
    assert.equal(empty.path.length,0);
    assert.equal(empty.turningPoints.length,0);
  }
  assert.equal(sim.configure({...base,casimir:.1,energy:1.5}).path.length,0);
});

test('Negative discriminant: no path, no curve, no current marker or fake diagnostics', () => {
  const sim=loadSimulator(), m=sim.configure({...base,casimir:.8,energy:1.40});
  near(m.roots.delta,-.1719);
  assert.equal(m.regime,'no real trajectory');
  assert.equal(m.path.length,0);
  assert.equal(m.radiusSeries.length,0);
  assert.equal(m.turningPoints.length,0);
  assert.equal(sim.run('energyCurvePoints(model).length'),0);
  assert.equal(m.canAnimate,false);
  assert.ok(sim.element('togglePlay').disabled);
  assert.equal(sim.element('trajectoryStatus').dataset.state,'empty');
  for (const id of ['metricPoint','metricCCheck','metricHCheck','metricResiduals']) {
    assert.match(sim.element(id).textContent,/no physical trajectory/);
  }
  sim.run('drawAll()');
  assert.equal(sim.drawCalls.filter(call=>call.method==='fill'&&call.color==='#f59e0b').length,0);
});

test('Animation freezes in empty and stationary modes and preserves pause intent', () => {
  const sim=loadSimulator();
  sim.run('animate(1000); animate(1050)');
  assert.ok(sim.run('currentPhase')>0);
  sim.element('togglePlay').dispatch('click');
  const paused=sim.run('currentPhase');
  sim.run('animate(1100)');
  assert.equal(sim.run('currentPhase'),paused);
  sim.element('preset-critical').dispatch('click');
  sim.run('animate(1150)');
  assert.equal(sim.run('currentPhase'),0);
  sim.configure({...base,casimir:.8,energy:1.4});
  sim.run('animate(1200)');
  assert.equal(sim.run('currentPhase'),0);
  sim.element('preset-guided').dispatch('click');
  assert.equal(sim.element('togglePlay').textContent,'play');
  sim.run('animate(1250)');
  assert.equal(sim.run('currentPhase'),0);
  sim.element('togglePlay').dispatch('click');
  sim.element('preset-critical').dispatch('click');
  sim.element('preset-singular').dispatch('click');
  assert.equal(sim.element('togglePlay').textContent,'pause');
  sim.run('animate(1300)');
  assert.ok(sim.run('currentPhase')>0);
  sim.element('resetParams').dispatch('click');
  assert.equal(sim.model().regime,'guided');
  assert.equal(sim.run('currentPhase'),0);
});

test('Keyboard controls: arrows, PageUp/Down, Home/End and preset state', () => {
  const sim=loadSimulator(), knob=sim.element('knob-casimir');
  knob.dispatch('keydown',{key:'Home'});
  assert.equal(sim.model().C,0);
  assert.equal(sim.model().regime,'singular cone');
  knob.dispatch('keydown',{key:'ArrowRight'});
  near(sim.model().C,.001);
  knob.dispatch('keydown',{key:'PageUp'});
  near(sim.model().C,.011);
  knob.dispatch('keydown',{key:'PageDown'});
  near(sim.model().C,.001);
  knob.dispatch('keydown',{key:'End'});
  assert.equal(sim.model().C,1);
  assert.equal(sim.model().path.length,0);
  assert.equal(sim.element('preset-guided').attributes['aria-pressed'],'false');
  assert.equal(knob.attributes['aria-valuenow'],'1.0000');
});

test('Guard rejects invalid states instead of clipping them into a trajectory', () => {
  const sim=loadSimulator();
  assert.equal(sim.run('isPhysicalState(model,[1,1,0])'),false);
  assert.equal(sim.run('isPhysicalState(model,[NaN,0,0])'),false);
  assert.equal(sim.run('isPhysicalState(model,[-1,0,0])'),false);
  sim.run('stateAtAngle=()=>[1,1,0]; refreshModel(true)');
  assert.equal(sim.model().path.length,0);
  assert.equal(sim.model().canAnimate,false);
  assert.match(sim.model().reason,/validation failed/);
});

test('Parameter grid covers control boundaries and independently predicts feasibility', () => {
  const sim=loadSimulator();
  let checked=0;
  for(const n0 of [1.05,1.5,1.8]) for(const alpha of [.01,.08,.22])
    for(const energy of [.4,1.3,1.4,1.5,1.58]) for(const casimir of [0,.001,.030625,.8,.98,1]) {
      // Avoid repeatedly drawing surfaces during this numeric sweep.
      sim.run(`Object.assign(controls.n0,{value:${n0}}); Object.assign(controls.alpha,{value:${alpha}}); Object.assign(controls.energy,{value:${energy}}); Object.assign(controls.casimir,{value:${casimir}}); model=computeModel();`);
      const m=sim.model(), b=n0*n0-energy*energy;
      const feasible=b>=0 && b*b>=4*alpha*casimir-1e-14;
      assert.equal(m.path.length>0,feasible,JSON.stringify({n0,alpha,energy,casimir}));
      if(feasible) allGeometry(sim,m);
      else assert.equal(sim.run('energyCurvePoints(model).length'),0);
      checked++;
    }
  assert.equal(checked,270);
});
