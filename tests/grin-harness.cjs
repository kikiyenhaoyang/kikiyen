const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

// Execute the production inline script, including its startup and event handlers.
// Canvas/DOM doubles are for logic tests only; they do not verify browser layout.
function loadSimulator(filename = path.join(__dirname, '../Web/poisson_variety_grin.html')) {
  const html = fs.readFileSync(filename, 'utf8');
  const scripts = [...html.matchAll(/<script\b[^>]*>([\s\S]*?)<\/script>/g)]
    .map(match => match[1]).filter(source => source.trim());
  const elements = new Map();
  const drawCalls = [];
  const context2D = new Proxy({
    createLinearGradient: () => ({addColorStop() {}}),
    measureText: text => ({width: text.length * 8}),
  }, {
    get(target, key) {
      if (key in target) return target[key];
      return (...args) => drawCalls.push({method: key, args, color: target.fillStyle});
    },
  });
  function element(id) {
    if (!elements.has(id)) {
      const handlers = new Map(), classes = new Set();
      elements.set(id, {
        id, width: 760, height: 760, textContent: '', innerHTML: '',
        disabled: false, dataset: {}, attributes: {}, style: {setProperty() {}},
        classList: {
          add: name => classes.add(name), remove: name => classes.delete(name),
          toggle: (name, force) => force ? classes.add(name) : classes.delete(name),
          contains: name => classes.has(name),
        },
        getContext: () => context2D,
        setAttribute(key, value) { this.attributes[key] = String(value); },
        addEventListener(name, fn) { handlers.set(name, fn); },
        setPointerCapture() {},
        getBoundingClientRect: () => ({left: 0, top: 0, width: 104, height: 104}),
        dispatch(name, values = {}) {
          if (name === 'click' && this.disabled) return;
          const event = {preventDefault() {}, ...values};
          handlers.get(name)?.(event);
        },
      });
    }
    return elements.get(id);
  }
  for (const name of ['guided', 'critical', 'singular']) element('preset-' + name);
  const sandbox = vm.createContext({
    console, document: {
      getElementById: element,
      querySelectorAll: () => [...elements.values()].filter(el => el.id.startsWith('preset-')),
    },
    requestAnimationFrame() {},
  });
  for (const script of scripts) vm.runInContext(script, sandbox, {filename});
  drawCalls.length = 0;
  const run = source => vm.runInContext(source, sandbox);
  return {
    run, element, drawCalls,
    model: () => run('model'),
    configure(values) {
      // Exact numeric inputs reproduce bug reports independently of knob rounding.
      sandbox.testValues = values;
      run('for(const [key,value] of Object.entries(testValues)) controls[key].value=value; refreshModel(true);');
      drawCalls.length = 0;
      return run('model');
    },
  };
}

function errors(model, points = model.path) {
  return {
    points: points.length,
    casimir: points.length ? Math.max(...points.map(([X,Y,Z]) => Math.abs(X*Y-Z*Z-model.C))) : null,
    energy: points.length ? Math.max(...points.map(([X,Y]) => Math.abs(model.n0*model.n0-model.alpha*X-Y-model.H*model.H))) : null,
  };
}
module.exports = {loadSimulator, errors};
