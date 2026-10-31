const {loadSimulator,errors} = require('./grin-harness.cjs');
if(!process.argv[2]) {
  console.error('Usage: node tests/grin-compare.cjs /path/to/baseline.html');
  process.exit(1);
}
const before=loadSimulator(process.argv[2]), after=loadSimulator();
const base={alpha:.08,n0:1.5,phase:0};
const cases=[
  ['Guided',{casimir:.030625,energy:1.4517}],
  ['Critical',{casimir:.98,energy:1.3}],
  ['Singular',{casimir:0,energy:1.42}],
  ['No real intersection',{casimir:.8,energy:1.4}],
  ['Singular with H>n0',{casimir:0,energy:1.58}],
];
console.log(JSON.stringify(cases.map(([name,values])=>{
  const parameters={...base,...values};
  const oldModel=before.configure(parameters), newModel=after.configure(parameters);
  return {name,parameters,before:errors(oldModel),after:errors(newModel),
    regime:newModel.regime,delta:newModel.roots.delta};
}),null,2));
