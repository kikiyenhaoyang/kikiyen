const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

// Execute the actual standalone page script with a minimal DOM/canvas adapter.
// This checks application logic and draw work, not browser rasterization/layout.
function loadPage(filename = path.join(__dirname, '../Web/satellite.html')) {
  const html = fs.readFileSync(filename, 'utf8');
  const calls = { lineTo: 0, moveTo: 0, alerts: [] };
  const ctx = new Proxy({
    measureText: text => ({ width: text.length * 8 }),
    createRadialGradient: () => ({ addColorStop() {} }),
  }, { get: (target, key) => key in target ? target[key] : (...args) => {
    if (key === 'lineTo' || key === 'moveTo') calls[key]++;
  }});
  const elements = new Map();
  for (const tag of html.matchAll(/<[^>]+\bid="([^"]+)"[^>]*>/g)) {
    const attr = name => tag[0].match(new RegExp(`\\b${name}="([^"]*)"`))?.[1];
    elements.set(tag[1], {
      value: attr('value') ?? '', textContent: '', hidden: /\bhidden\b/.test(tag[0]),
      width: Number(attr('width')), height: Number(attr('height')),
      attributes: {}, listeners: {}, validity: { badInput: false },
      get valueAsNumber() { return this.value.trim() === '' ? NaN : Number(this.value); },
      getContext: () => ctx,
      addEventListener(event, fn) { this.listeners[event] = fn; },
      setAttribute(key, value) { this.attributes[key] = value; },
      removeAttribute(key) { delete this.attributes[key]; },
      setCustomValidity(message) { this.validationMessage = message; },
      focus() { this.focused = true; },
    });
  }
  const sandbox = {
    document: { getElementById: id => elements.get(id) },
    requestAnimationFrame() {}, alert: message => calls.alerts.push(message),
  };
  const context = vm.createContext(sandbox);
  for (const script of html.matchAll(/<script\b[^>]*>([\s\S]*?)<\/script>/g)) {
    if (script[1].trim()) vm.runInContext(script[1], context, { filename, timeout: 10000 });
  }
  return {
    html, elements, calls, context,
    run: expression => vm.runInContext(expression, context, { timeout: 10000 }),
    params: () => sandbox.readParams(),
    simulate: params => sandbox.simulate(params),
    draw: (k, data) => sandbox.drawAll(k, data),
  };
}

module.exports = { loadPage };
