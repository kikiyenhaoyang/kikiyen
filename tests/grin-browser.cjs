const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const path = require('node:path');
const http = require('node:http');
const {chromium} = require('playwright');

// Real Chromium tests served over loopback HTTP; no DOM/canvas doubles.
// Mobile is viewport/touch emulation, not a physical device or Mobile Safari.
async function main() {
  const html = await fs.readFile(path.join(__dirname,'../Web/poisson_variety_grin.html'));
  const artifacts = process.env.GRIN_ARTIFACT_DIR || path.join(__dirname,'artifacts');
  await fs.mkdir(artifacts,{recursive:true});
  const server=http.createServer((req,res)=>{
    if(req.url!=='/'){res.writeHead(404);res.end();return;}
    res.writeHead(200,{'Content-Type':'text/html; charset=utf-8'}); res.end(html);
  });
  await new Promise(resolve=>server.listen(0,'127.0.0.1',resolve));
  let browser;
  const report={engine:'Chromium',results:[]};
  try {
    browser=await chromium.launch();
    report.version=browser.version();
    for(const device of [
      {name:'desktop',viewport:{width:1440,height:1000},isMobile:false,hasTouch:false},
      {name:'mobile',viewport:{width:390,height:844},isMobile:true,hasTouch:true},
    ]) {
      const {name,...options}=device;
      const context=await browser.newContext({...options,deviceScaleFactor:1});
      const page=await context.newPage(), pageErrors=[];
      page.on('pageerror',error=>pageErrors.push(error.message));
      await page.goto(`http://127.0.0.1:${server.address().port}/`,{waitUntil:'networkidle'});
      const state=()=>page.evaluate(()=>({
        regime:model.regime,points:model.path.length,canAnimate:model.canAnimate,
        C:model.C,H:model.H,phase:currentPhase,playing,
        maxC:model.maxCError,maxH:model.maxHError,roots:model.turningPoints,
        path:model.path,series:model.radiusSeries.map(({r})=>r),
      }));
      const valid=async()=>{
        const s=await state();
        assert.ok(s.points>0);
        assert.ok(s.maxC<1e-12 && s.maxH<1e-12);
        assert.ok(s.path.every(([X,Y,Z])=>X>=0&&Y>=0&&[X,Y,Z].every(Number.isFinite)));
        return s;
      };
      const shot=async suffix=>{
        await page.locator('#simulator > .card').screenshot({path:path.join(artifacts,`${name}-${suffix}.png`)});
      };
      const layout=await page.evaluate(expectedWidth=>({
        width:innerWidth,scroll:document.documentElement.scrollWidth,
        katexLoaded:typeof katex!=='undefined',
        overflow:[...document.querySelectorAll('*')].filter(el=>{
          const r=el.getBoundingClientRect();
          if(r.width<=1 || r.right<=expectedWidth+1 || el.closest('.katex-mathml')) return false;
          for(let parent=el.parentElement;parent;parent=parent.parentElement){
            if(['hidden','clip','auto','scroll'].includes(getComputedStyle(parent).overflowX)
              && parent.getBoundingClientRect().right<=expectedWidth+1) return false;
          }
          return true;
        }).slice(0,30).map(el=>({
          tag:el.tagName,class:el.className,right:el.getBoundingClientRect().right,
          text:el.textContent.slice(0,120),
          parent:el.parentElement?.className,
          parentOverflow:el.parentElement?getComputedStyle(el.parentElement).overflowX:null,
        })),
      }),device.viewport.width);
      await page.screenshot({path:path.join(artifacts,name+'-page.png'),fullPage:true});
      // Mobile browsers can expand innerWidth to fit overflowing content.
      // Comparing scrollWidth only with innerWidth would falsely pass that case.
      assert.equal(layout.width,device.viewport.width,JSON.stringify(layout));
      assert.ok(layout.scroll<=device.viewport.width+1,JSON.stringify(layout));
      assert.equal((await valid()).regime,'guided');
      await page.getByRole('button',{name:'pause',exact:true}).click();
      const paused=(await state()).phase;
      await page.waitForTimeout(100);
      assert.equal((await state()).phase,paused);
      await shot('guided');

      // Navigate between presets with the real keyboard, not a click handler stub.
      await page.locator('#preset-guided').focus();
      await page.keyboard.press('Tab');
      assert.equal(await page.evaluate(()=>document.activeElement.id),'preset-critical');
      await page.keyboard.press('Enter');
      let s=await valid();
      assert.equal(s.regime,'critical'); assert.equal(s.points,1);
      assert.equal(s.canAnimate,false);
      assert.equal(new Set(s.series).size,1);
      assert.equal(await page.locator('#togglePlay').isDisabled(),true);
      assert.match(await page.locator('#trajectoryStatus').innerText(),/fixed reduced state/);
      await shot('critical');

      // Reproduce the reported C=.8, H=1.40 case through keyboard controls.
      await page.locator('#knob-casimir').focus();
      await page.keyboard.press('End');
      for(let i=0;i<20;i++) await page.keyboard.press('PageDown');
      await page.locator('#knob-energy').focus();
      for(let i=0;i<10;i++) await page.keyboard.press('PageUp');
      s=await state();
      assert.ok(Math.abs(s.C-.8)<1e-12 && Math.abs(s.H-1.4)<1e-12);
      assert.equal(s.points,0); assert.equal(s.phase,0);
      assert.equal(await page.locator('#togglePlay').isDisabled(),true);
      assert.match(await page.locator('#metricPoint').innerText(),/no physical trajectory/);
      assert.match(await page.locator('#metricHCheck').innerText(),/no physical trajectory/);
      assert.equal(await page.evaluate(()=>energyCurvePoints(model).length),0);
      await shot('empty');

      await page.getByRole('button',{name:'guided',exact:true}).click();
      assert.equal(await page.locator('#togglePlay').innerText(),'play');
      await page.getByRole('button',{name:'play',exact:true}).click();
      await page.waitForFunction(()=>currentPhase>0);
      await page.getByRole('button',{name:'critical',exact:true}).click();
      await page.waitForTimeout(100);
      assert.equal((await state()).phase,0);
      await page.getByRole('button',{name:'C=0 cone',exact:true}).click();
      await page.waitForFunction(()=>currentPhase>0);
      s=await valid();
      assert.equal(s.regime,'singular cone'); assert.equal(s.roots[0],0);
      await shot('singular');
      await page.getByRole('button',{name:'pause',exact:true}).click();

      await page.locator('#knob-casimir').focus();
      await page.keyboard.press('Home'); await page.keyboard.press('ArrowRight');
      assert.equal(await page.locator('#knob-casimir').getAttribute('aria-valuenow'),'0.0010');
      assert.equal(await page.locator('#preset-guided').getAttribute('aria-pressed'),'false');
      await page.keyboard.press('ArrowLeft');
      assert.equal((await state()).C,0);
      await page.keyboard.press('End');
      assert.equal((await state()).points,0);

      // Exercise pointer controls with mouse on desktop and touch on mobile.
      await page.locator('#knob-casimir').scrollIntoViewIfNeeded();
      const knob=await page.locator('#knob-casimir').boundingBox();
      if(device.hasTouch) await page.touchscreen.tap(knob.x+knob.width/2,knob.y+5);
      else await page.mouse.click(knob.x+knob.width/2,knob.y+5);
      assert.ok(Math.abs((await state()).C-.5)<.002);
      await page.getByRole('button',{name:'guided defaults',exact:true}).click();
      await valid();
      if(!device.hasTouch){
        await page.locator('#varietyCanvas').scrollIntoViewIfNeeded();
        const canvas=await page.locator('#varietyCanvas').boundingBox();
        const yaw=await page.evaluate(()=>cameraYaw);
        await page.mouse.move(canvas.x+canvas.width/2,canvas.y+canvas.height/2);
        await page.mouse.down();
        await page.mouse.move(canvas.x+canvas.width/2+40,canvas.y+canvas.height/2+20,{steps:5});
        await page.mouse.up();
        assert.notEqual(await page.evaluate(()=>cameraYaw),yaw);
        await page.getByRole('button',{name:'reset view',exact:true}).click();
        assert.ok(Math.abs(await page.evaluate(()=>cameraYaw)-yaw)<1e-12);
      }
      assert.deepEqual(pageErrors,[]);
      report.results.push({device:name,viewport:device.viewport,layout,status:'passed',
        checks:['guided','critical','empty','singular','pause/resume','Tab/Enter',
          'arrows/PageUp/PageDown/Home/End',device.hasTouch?'touch knob':'mouse knob/camera',
          'no horizontal overflow','no uncaught page errors']});
      await context.close();
    }
    console.log(JSON.stringify(report,null,2));
  } finally {
    await fs.writeFile(path.join(artifacts,'browser-results.json'),JSON.stringify(report,null,2));
    if(browser) await browser.close();
    await new Promise(resolve=>server.close(resolve));
  }
}
main().catch(error=>{console.error(error);process.exitCode=1;});
