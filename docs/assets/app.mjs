import {sensitivityParameters} from './sensitivity.mjs?v=sensitivity-1';
import {D,mod,inverseMoll,moll,vectorPixel,fitsMap,fromEquatorial,setUTCField,shiftUTCField,ssoElements} from './sky-core.mjs?v=sensitivity-1';
const isAir=kind=>['airglow','total','sensitivity'].includes(kind);
const $=id=>document.getElementById(id),canvas=$('map'),ctx=canvas.getContext('2d');
let manifest,worker,result,savedRequest,date=new Date(),generation=0;
const fields=['year','month','day','hour','minute','second'],labels=['년','월','일','시','분','초'];
const getters=['getUTCFullYear','getUTCMonth','getUTCDate','getUTCHours','getUTCMinutes','getUTCSeconds'];
function dateUI(){fields.forEach((f,i)=>{$(f).value=String(date[getters[i]]()+(i===1?1:0)).padStart(i===0?4:2,'0');});}
fields.forEach((f,i)=>{const wrap=document.createElement('div');wrap.className='datepart';wrap.innerHTML=`<input id="${f}" type="text" inputmode="numeric" aria-label="UTC ${labels[i]}"><div class="date-arrows"><button type="button" aria-label="${labels[i]} 증가">▲</button><button type="button" aria-label="${labels[i]} 감소">▼</button></div><span>${labels[i]}</span>`;$('datetime').append(wrap);
 function change(delta){date=shiftUTCField(date,i,delta);dateUI();schedule();}
 wrap.querySelectorAll('button').forEach((b,j)=>b.onclick=()=>change(j===0?1:-1));$(f).onkeydown=e=>{if(e.key==='ArrowUp'||e.key==='ArrowDown'){e.preventDefault();change(e.key==='ArrowUp'?1:-1);}};
 $(f).onchange=()=>{const v=Number($(f).value);if(Number.isInteger(v)){date=setUTCField(date,i,v-(i===1?1:0));}dateUI();schedule();};});
date.setUTCMilliseconds(0);dateUI();$('now').onclick=()=>{date=new Date();date.setUTCMilliseconds(0);dateUI();schedule();};
let timer;function schedule(){clearTimeout(timer);timer=setTimeout(calculate,250);}
function panels(){const sensitivity=$('kind').value==='sensitivity';$('sensitivity-controls').hidden=!sensitivity;
 for(const id of ['anti','moon-avoid']){if(sensitivity)$(id).checked=true;$(id).disabled=sensitivity;}
 if(sensitivity)$('model').value='combined';$('model').disabled=sensitivity;
 $('moon-controls').hidden=false;$('zodi-controls').hidden=false;$('epoch').disabled=!$('evolve').checked;$('dgl-hint').hidden=!['dgl','zodiacal_dgl','total','sensitivity'].includes($('kind').value);$('air-controls').hidden=false;$('air-help').hidden=false;$('air-model-control').hidden=!isAir($('kind').value);$('air-component-help').hidden=!isAir($('kind').value);}
for(const id of ['kind','band','coords','nside','anti','elongation','moon-avoid','moon-angle','height','orbit','phase','min-elv','evolve','epoch','model','f-sys','exposure-hours','bin-side','snr','n-reads'])$(id).onchange=()=>{panels();schedule();};
// Sun-synchronous inclination for the entered height (same J2 formula as the orbit model).
function showInclination(){const H=Number($('height').value);let text='—';try{if(H>0&&H<1500)text=ssoElements(H).inclinationDeg.toFixed(2);}catch{}$('inclination').textContent=text;}
$('height').addEventListener('input',showInclination);showInclination();
for(const id of ['contour','contour-step'])$(id).onchange=()=>{if(result)draw();};$('calculate').onclick=calculate;
// One long-lived worker keeps its decompressed, SHA-checked arrays; each request carries an id and
// the worker abandons a request as soon as a newer id or a cancel message arrives.
function idle(){$('cancel').disabled=true;$('progress').hidden=true;}
function stop(){generation++;worker?.postMessage({cancel:true});idle();}
function ensureWorker(){if(worker)return worker;worker=new Worker(new URL('./worker.mjs?v=sensitivity-1',import.meta.url),{type:'module'});worker.onerror=e=>{worker?.terminate();worker=null;generation++;idle();$('status').textContent='계산 실패: '+e.message;};return worker;}
$('cancel').onclick=()=>{clearTimeout(timer);stop();$('status').textContent='계산을 중지했습니다. 표시 중인 지도는 이전 결과입니다.';};
function request(){const r={kind:$('kind').value,coords:$('coords').value,band:$('band').value,nside:Number($('nside').value),utc:date.toISOString(),anti:$('anti').checked,elongation:Number($('elongation').value),moonAvoid:$('moon-avoid').checked,moonAngle:Number($('moon-angle').value),height:Number($('height').value),orbit:$('orbit').value,orbitalPhaseDeg:Number($('phase').value),minElvDeg:Number($('min-elv').value),evolve:$('evolve').checked,epoch:$('epoch').value,model:$('model').value};
 if(r.kind==='sensitivity')r.sensitivity=sensitivityParameters(r.band,{f_sys:Number($('f-sys').value)/100,exposure_hours:Number($('exposure-hours').value),bin_side:Number($('bin-side').value),snr:Number($('snr').value),n_reads:Number($('n-reads').value)});
 if(!Number.isFinite(r.moonAngle)||r.moonAngle<0||r.moonAngle>180)throw Error('Moon separation은 0–180°입니다.');if(r.elongation<0||r.elongation>180||!Number.isFinite(r.elongation))throw Error('Sun separation은 0–180°입니다.');if((!(r.height>0&&r.height<1500)||!Number.isFinite(r.orbitalPhaseDeg)))throw Error('H는 0–1500 km 사이이며 위상은 유한한 값이어야 합니다.');if((!Number.isFinite(r.minElvDeg)||r.minElvDeg < -90||r.minElvDeg>180))throw Error('최소 ELV는 −90–180°입니다.');if(r.evolve){if(!/(Z|[+-]\d\d:\d\d)$/.test(r.epoch))r.epoch+='Z';if(!Number.isFinite(Date.parse(r.epoch)))throw Error('위상 기준 UTC를 확인하십시오.');}return r;}
function calculate(){clearTimeout(timer);if(!manifest)return;try{const r=request();stop();const id=generation;ensureWorker();$('cancel').disabled=false;$('progress').hidden=false;$('progress').value=0;$('status').textContent='데이터 읽기 및 계산 중…';const start=performance.now();
 worker.onmessage=({data:d})=>{if(d.id!==id||id!==generation)return;if(d.error){idle();$('status').textContent=d.error;return;}if(d.progress!==undefined){$('progress').value=d.progress;$('status').textContent=`계산 중 · ${Math.round(d.progress*100)}%`;return;}result=d.result;savedRequest=r;idle();draw();$('status').textContent=`완료 · ${result.flux.length.toLocaleString()} pixels · ${((performance.now()-start)/1000).toFixed(1)} s${result.view&&!result.view.some(Boolean)?' · 관측 영역 없음: 위상·ELV·회피 조건을 확인하십시오.':''}`;for(const k of ['fits','pdf','png'])$(k).disabled=false;};
 worker.postMessage({id,request:r,manifest});}catch(e){$('status').textContent=e.message;}}
const area={x:90,y:52,w:1020,h:510},gridW=340,gridH=170;
function sample(v){const i=vectorPixel(savedRequest.nside,v);return {mag:result.mag[i],i};}
const stops=[[253,231,37],[94,201,98],[33,145,140],[59,82,139],[68,1,84]];
function color(v,lo,hi){const u=Math.max(0,Math.min(4,4*(v-lo)/(hi-lo))),i=Math.min(3,Math.floor(u)),f=u-i;return stops[i].map((x,j)=>Math.round(x*(1-f)+stops[i+1][j]*f));}
let scale;
function coordinateLabels(){return savedRequest.coords==='galactic'?'l, b':'RA, DEC';}
function projectSky(v){const p=fromEquatorial(v,savedRequest.coords);return moll(Math.atan2(p[1],p[0]),Math.asin(Math.max(-1,Math.min(1,p[2]))));}
function bodyCoordinates(v){const p=fromEquatorial(v,savedRequest.coords);return `(${mod(Math.atan2(p[1],p[0])/D,360).toFixed(2)}°, ${(Math.asin(Math.max(-1,Math.min(1,p[2])))/D).toFixed(2)}°)`;}
function draw(){$('cursor').textContent=`— mag at (${coordinateLabels()}) = (—, —)`;canvas.setAttribute('aria-label',coordinateLabels()+' Mollweide sky map');const r=savedRequest,b=manifest.bands[r.band];$('map-title').textContent=`${r.coords==='galactic'?'Galactic':'Equatorial'} | ${({zodiacal:'Zodiacal light',dgl:'DGL',zodiacal_dgl:'Zodiacal light + DGL',airglow:'Airglow',total:'Zodiacal light + DGL + Airglow',sensitivity:'Sensitivity (Zodiacal + DGL + Airglow)'})[r.kind]} | ${r.band} · ${b[0]}–${b[1]} nm | ${'UTC '+r.utc.replace('T',' ').slice(0,19)+' · '+r.orbit+' · H = '+r.height+' km · φ = '+result.state.orbitalPhaseDeg.toFixed(1)+'° · SZA = '+result.state.sza.toFixed(1)+'° · ELV ≥ '+r.minElvDeg+'°'} | NSIDE ${r.nside}${result.moon&&r.moonAvoid?' · Moon separation ≥ '+r.moonAngle+'°':''}${r.anti?' · Sun separation ≥ '+r.elongation+'°':''}`;
 if(r.kind==='sensitivity'){const p=r.sensitivity;$('map-title').textContent+=` | ${p.snr}σ · ${p.exposure_hours} h usable · f_sys=${100*p.f_sys}% · ${p.bin_side}×${p.bin_side} pixels · ${p.n_reads} reads`;}
 const finite=Array.from(result.mag).filter(Number.isFinite);let lo=Infinity,hi=-Infinity;for(const x of finite){lo=Math.min(lo,x);hi=Math.max(hi,x);}if(!finite.length){lo=18;hi=24;}lo=Math.floor(lo*2)/2;hi=Math.max(lo+.5,Math.ceil(hi*2)/2);scale={lo,hi};
 ctx.fillStyle='white';ctx.fillRect(0,0,canvas.width,canvas.height);const image=ctx.createImageData(area.w,area.h);
 for(let y=0;y<area.h;y++)for(let x=0;x<area.w;x++){const v=inverseMoll(2*(x+.5)/area.w-1,1-2*(y+.5)/area.h);if(!v)continue;const {mag,i}=sample(v);let rgb=Number.isFinite(mag)?color(mag,lo,hi):[245,245,245];if(result.earth?.[i]&&(!result.view||result.view[i]))rgb=[174,181,190];else if(result.valid&&(!result.view||result.view[i])&&!result.valid[i]&&(x+y)%12<2)rgb=[110,110,110];const p=(y*area.w+x)*4;image.data.set([...rgb,255],p);}ctx.putImageData(image,area.x,area.y);
 ctx.strokeStyle='#ffffff60';ctx.lineWidth=1;
 function path(points){ctx.beginPath();let first=true;for(const [x,y]of points){const px=area.x+(x+1)*area.w/2,py=area.y+(1-y)*area.h/2;if(first){ctx.moveTo(px,py);first=false;}else ctx.lineTo(px,py);}ctx.stroke();}
 for(const dec of [-60,-30,0,30,60])path(Array.from({length:181},(_,i)=>moll((-180+i*2)*D,dec*D)));
 for(const ra of [-150,-120,-90,-60,-30,0,30,60,90,120,150])path(Array.from({length:181},(_,i)=>moll(ra*D,(-90+i)*D)));
 if($('contour').checked)contours(lo,hi);
 drawMoon();drawSun();drawZenith();ctx.lineWidth=1;ctx.strokeStyle='#5b6b7a';ctx.beginPath();ctx.ellipse(area.x+area.w/2,area.y+area.h/2,area.w/2,area.h/2,0,0,2*Math.PI);ctx.stroke();
 ctx.fillStyle='#203348';ctx.font='14px system-ui';ctx.textAlign='center';for(const ra of [150,90,30,0,330,270,210]){const [x]=moll(ra*D,0);ctx.fillText(ra+'°',area.x+(x+1)*area.w/2,area.y+area.h+24);}ctx.fillText(r.coords==='galactic'?'l (Galactic longitude) — increases to the left':'RA (Equatorial) — increases to the left',600,area.y+area.h+48);ctx.save();ctx.translate(25,area.y+area.h/2);ctx.rotate(-Math.PI/2);ctx.fillText(r.coords==='galactic'?'b (Galactic latitude, °)':'DEC (°)',0,0);ctx.restore();for(const dec of [-60,-30,0,30,60]){const [,y]=moll(0,dec*D);ctx.textAlign='right';ctx.fillText(dec+'°',area.x-12,area.y+(1-y)*area.h/2+5);}
 for(let x=0;x<700;x++){ctx.fillStyle=`rgb(${color(lo+(hi-lo)*x/699,lo,hi)})`;ctx.fillRect(250+x,644,1,18);}ctx.fillStyle='#203348';ctx.textAlign='center';for(let i=0;i<=6;i++){const x=250+700*i/6;ctx.fillText((lo+(hi-lo)*i/6).toFixed(1),x,681);}ctx.fillText(r.kind==='sensitivity'?'Limiting AB mag / arcsec² (larger = deeper)':'AB mag / arcsec²',600,707);
 $('legend').textContent=`Contour: ${Number($('contour-step').value).toFixed(1)} mag · ${lo.toFixed(1)}–${hi.toFixed(1)} mag / arcsec²${result.moon?' · Moon (from satellite, '+coordinateLabels()+') = '+bodyCoordinates(result.moon):''}`;}
function drawSeparation(v,angle,color){
 if(!(angle>0&&angle<180))return;
 const project=p=>{const q=projectSky(p);return [area.x+(q[0]+1)*area.w/2,area.y+(1-q[1])*area.h/2];};
 ctx.save();ctx.strokeStyle=color;ctx.lineWidth=1.4;
  const reference=Math.abs(v[2])<.9?[0,0,1]:[1,0,0];let u=[v[1]*reference[2]-v[2]*reference[1],v[2]*reference[0]-v[0]*reference[2],v[0]*reference[1]-v[1]*reference[0]];const n=Math.hypot(...u);u=u.map(x=>x/n);
  const w=[v[1]*u[2]-v[2]*u[1],v[2]*u[0]-v[0]*u[2],v[0]*u[1]-v[1]*u[0]],a=angle*D;let previous;
  ctx.setLineDash([5,4]);ctx.beginPath();for(let i=0;i<=360;i++){const t=i*D,p=project(v.map((x,j)=>Math.cos(a)*x+Math.sin(a)*(Math.cos(t)*u[j]+Math.sin(t)*w[j])));if(!previous||Math.hypot(p[0]-previous[0],p[1]-previous[1])>area.w/4)ctx.moveTo(...p);else ctx.lineTo(...p);previous=p;}ctx.stroke();ctx.setLineDash([]);
 ctx.restore();
}
function drawZenith(){
 if(!result.state)return;const v=result.state.zenith;
 drawSeparation(v,result.state.limbZenithDeg-savedRequest.minElvDeg,'#008060');
 const q=projectSky(v),x=area.x+(q[0]+1)*area.w/2,y=area.y+(1-q[1])*area.h/2;
 ctx.save();ctx.fillStyle='#006bbb';ctx.beginPath();ctx.arc(x,y,4,0,Math.PI*2);ctx.fill();ctx.font='bold 12px system-ui';ctx.textAlign=x>area.x+area.w-65?'right':'left';const dx=ctx.textAlign==='right'?-7:7;ctx.strokeStyle='white';ctx.lineWidth=3;ctx.strokeText('Zenith',x+dx,y+15);ctx.fillText('Zenith',x+dx,y+15);ctx.restore();
}
function drawSun(){
 if(!result.sun)return;
 if(savedRequest.anti)drawSeparation(result.sun,savedRequest.elongation,'#e83131');
 const v=result.sun,q=projectSky(v),x=area.x+(q[0]+1)*area.w/2,y=area.y+(1-q[1])*area.h/2;
 ctx.save();ctx.beginPath();for(let i=0;i<10;i++){const a=-Math.PI/2+i*Math.PI/5,r=i%2?3.5:8,px=x+r*Math.cos(a),py=y+r*Math.sin(a);if(i)ctx.lineTo(px,py);else ctx.moveTo(px,py);}ctx.closePath();ctx.fillStyle='#e83131';ctx.fill();ctx.lineWidth=1;ctx.strokeStyle='#203348';ctx.stroke();
 ctx.font='bold 12px system-ui';ctx.textAlign=x>area.x+area.w-55?'right':'left';ctx.textBaseline='bottom';const lx=ctx.textAlign==='right'?x-10:x+10;ctx.lineWidth=3;ctx.strokeStyle='white';ctx.strokeText('Sun',lx,y-9);ctx.fillStyle='#c02020';ctx.fillText('Sun',lx,y-9);ctx.restore();
}
function drawMoon(){
 if(!result.moon)return;const v=result.moon,r=savedRequest;
 const project=p=>{const q=projectSky(p);return [area.x+(q[0]+1)*area.w/2,area.y+(1-q[1])*area.h/2];};
 ctx.save();ctx.strokeStyle='#b56d00';ctx.lineWidth=1.4;
 if(r.moonAvoid)drawSeparation(v,r.moonAngle,'#b56d00');
 const [x,y]=project(v);ctx.beginPath();ctx.arc(x,y,6,0,2*Math.PI);ctx.fillStyle='#ffcf56';ctx.fill();ctx.strokeStyle='#203348';ctx.stroke();ctx.font='bold 12px system-ui';ctx.textAlign=x>area.x+area.w-65?'right':'left';ctx.textBaseline='middle';const lx=ctx.textAlign==='right'?x-10:x+10;ctx.strokeStyle='white';ctx.lineWidth=3;ctx.strokeText('Moon',lx,y);ctx.fillStyle='#805000';ctx.fillText('Moon',lx,y);ctx.restore();
}
function contours(lo,hi){const w=gridW,h=gridH,a=new Float32Array((w+1)*(h+1));a.fill(NaN);for(let y=0;y<=h;y++)for(let x=0;x<=w;x++){const v=inverseMoll(2*x/w-1,1-2*y/h);if(v){const s=sample(v);if(!result.valid||result.valid[s.i])a[y*(w+1)+x]=s.mag;}}
 ctx.lineWidth=.8;ctx.strokeStyle='#ffffffb0';ctx.font='11px system-ui';const usedLabels=[];
 const step=Number($('contour-step').value),first=Math.floor((lo+1e-9)/step)+1,last=Math.ceil((hi-1e-9)/step)-1,count=Math.max(1,last-first+1);
 for(let k=first;k<=last;k++){const level=Number((k*step).toFixed(1));ctx.lineWidth=count>1?.6+1.6*(k-first)/(count-1):1.2;ctx.beginPath();let label=null;for(let y=0;y<h;y++)for(let x=0;x<w;x++){const vals=[a[y*(w+1)+x],a[y*(w+1)+x+1],a[(y+1)*(w+1)+x+1],a[(y+1)*(w+1)+x]];if(!vals.every(Number.isFinite))continue;const corners=[[x,y],[x+1,y],[x+1,y+1],[x,y+1]],p=[];
 for(let k=0;k<4;k++){const j=(k+1)%4;if((vals[k]<level)!==(vals[j]<level)){const f=(level-vals[k])/(vals[j]-vals[k]);p.push([area.x+(corners[k][0]+f*(corners[j][0]-corners[k][0]))*area.w/w,area.y+(corners[k][1]+f*(corners[j][1]-corners[k][1]))*area.h/h]);}}
 // Resolve saddle cells with their bilinear center value.
 if(p.length===4&&((vals.reduce((s,v)=>s+v,0)/4<level)!==(vals[0]<level)))p.push(p.shift());
 for(let k=0;k+1<p.length;k+=2){ctx.moveTo(...p[k]);ctx.lineTo(...p[k+1]);if(!label&&x>w*.15&&x<w*.85&&y>h*.15&&y<h*.85&&usedLabels.every(a=>Math.hypot(a[0]-p[k][0],a[1]-p[k][1])>30))label=p[k];}}
 ctx.stroke();if(label){usedLabels.push(label);ctx.fillStyle='#203348';ctx.fillRect(label[0]-13,label[1]-9,26,13);ctx.fillStyle='white';ctx.textAlign='center';ctx.fillText(level.toFixed(1),label[0],label[1]+1);}}
}
canvas.onmousemove=e=>{if(!result)return;const box=canvas.getBoundingClientRect(),x=(e.clientX-box.left)*canvas.width/box.width,y=(e.clientY-box.top)*canvas.height/box.height,v=inverseMoll(2*(x-area.x)/area.w-1,1-2*(y-area.y)/area.h);if(!v){$('cursor').textContent=`— mag at (${coordinateLabels()}) = (—, —)`;return;}const {mag,i}=sample(v),ra=mod(Math.atan2(v[1],v[0])/D,360),dec=Math.asin(v[2])/D;$('cursor').textContent=`${Number.isFinite(mag)?mag.toFixed(3):mag===Infinity?'∞':'—'}${savedRequest.kind==='sensitivity'?' limiting':''} mag at (${coordinateLabels()}) = (${ra.toFixed(2)}°, ${dec.toFixed(2)}°)${result.theta?' · θ = '+result.theta[i].toFixed(2)+'° · ELV = '+result.elv[i].toFixed(2)+'°':''}${result.view&&!result.view[i]?' · excluded by avoidance / ELV':result.earth?.[i]?' · Earth occulted':result.valid&&!result.valid[i]?' · outside model support':''}`;};
function filename(){const r=savedRequest,stamp=r.utc.slice(0,19).replaceAll('-','').replaceAll(':','').replace('T','_')+`_${r.orbit}_H${r.height}_phi${result.state.orbitalPhaseDeg.toFixed(3)}_elv${r.minElvDeg}${isAir(r.kind)?'_'+r.model:''}`;return `${r.kind}${r.kind==='sensitivity'?`_t${r.sensitivity.exposure_hours}h_fsys${r.sensitivity.f_sys}_bin${r.sensitivity.bin_side}_snr${r.sensitivity.snr}_reads${r.sensitivity.n_reads}`:''}_${r.band}_nside${r.nside}_${stamp}_${r.coords}${result.moon&&r.moonAvoid?'_moon'+r.moonAngle:''}${r.anti?'_anti-sun':''}`;}
function download(blob,name){const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),10000);}
$('fits').onclick=async()=>{try{const bytes=fitsMap(result,savedRequest,manifest),blob=await new Response(new Blob([bytes]).stream().pipeThrough(new CompressionStream('gzip'))).blob();download(blob,filename()+'.fits.gz');}catch(e){$('status').textContent='FITS 저장 실패: '+e.message;}};
$('pdf').onclick=()=>{document.title=filename();window.print();};$('png').onclick=()=>canvas.toBlob(blob=>download(blob,filename()+'.png'));
try{if(!('Worker'in window)||!('DecompressionStream'in window)||!('CompressionStream'in window))throw Error('최신 Chrome, Edge, Firefox 또는 Safari가 필요합니다.');const response=await fetch(new URL('../data/manifest.json?v=sensitivity-1',import.meta.url));if(!response.ok)throw Error('manifest 읽기 실패');manifest=await response.json();$('moon-avoid').checked=true;$('moon-angle').value='60';$('anti').checked=true;$('elongation').value='90';panels();calculate();}catch(e){$('status').textContent=e.message+' · HTML을 직접 열지 말고 HTTP 서버 또는 GitHub Pages에서 실행하십시오.';}
