// Optional real-browser smoke test: requires Playwright and its Chromium binary.
// Run: node tests/satellite.browser.cjs
const assert = require('node:assert/strict');
const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const { chromium } = require('playwright');

(async () => {
  const html = fs.readFileSync(path.join(__dirname, '../Web/satellite.html'));
  const server = http.createServer((req, res) => {
    res.writeHead(200, { 'Content-Type': 'text/html; charset=utf-8' });
    res.end(html);
  });
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  let browser;
  try {
    browser = await chromium.launch({ headless: true });
    for (const width of [1280, 390]) {
      const page = await browser.newPage({ viewport: { width, height: 900 } });
      const errors = [];
      page.on('pageerror', error => errors.push(error.message));
      // The external math typesetter is unrelated to the local simulator tests.
      await page.route('https://cdn.jsdelivr.net/**', route => route.abort());
      await page.goto(`http://127.0.0.1:${server.address().port}/satellite.html`);
      assert.equal(await page.evaluate(() => simData.nframes), 720);
      assert.equal(await page.locator('#metricVisible').textContent(), '100.0%');
      const input = page.locator('#samples');
      const restart = page.locator('#restartBtn');
      await page.evaluate(() => { window.previousData = simData; });
      for (const value of ['', 'NaN', 'Infinity', '1e309', '720.5', '119', '10001', '150000']) {
        // Assign invalid text via the DOM to exercise native number sanitization.
        await input.evaluate((el, next) => { el.value = next; }, value);
        await restart.click();
        assert.equal(await page.evaluate(() => simData === window.previousData), true);
        assert.equal(await input.getAttribute('aria-invalid'), 'true');
        assert.equal(await page.locator('#simulationError').isVisible(), true);
        assert.equal(await input.evaluate(el => el === document.activeElement), true);
      }
      for (const value of ['120', '121', '7.2e2', '10000']) {
        await input.fill(value);
        await restart.click();
        assert.equal(await page.evaluate(() => simData.nframes), Number(value));
        assert.equal(await input.evaluate(el => el.checkValidity()), true);
        assert.equal(await page.locator('#simulationError').isVisible(), false);
        assert.equal(await input.getAttribute('aria-invalid'), null);
      }
      await page.waitForFunction(() => frameIndex > 0);
      assert.equal(await page.evaluate(() => ['earthCanvas', 'groundCanvas', 'skyCanvas'].every(id => {
        const canvas = document.getElementById(id);
        return canvas.getContext('2d').getImageData(0, 0, canvas.width, canvas.height).data
          .some((channel, index) => index % 4 !== 3 && channel > 0);
      })), true);
      assert.deepEqual(errors, []);
      console.log(`PASS Chromium: ${width}px, defaults, invalid inputs, recovery, 10,000 samples and all canvases`);
      await page.close();
    }
  } finally {
    if (browser) await browser.close();
    await new Promise(resolve => server.close(resolve));
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
