import {contourPaths,contourLabels} from './contours.mjs?v=contour-labels-1';
import {D,mod,inverseMoll,moll,vectorPixel,fitsMap} from './sky-core.mjs?v=dgl-selfabs-1';
const $=id=>document.getElementById(id),canvas=$('map'),ctx=canvas.getContext('2d');
let manifest,worker,result,savedRequest,date=new Date(),generation=0;
const fields=['year','month','day','hour','minute','second'],labels=['년','월','일','시','분','초'];
const getters=['getUTCFullYear','getUTCMonth','getUTCDate','getUTCHours','getUTCMinutes','getUTCSeconds'],setters=['setUTCFullYear','setUTCMonth','setUTCDate','setUTCHours','setUTCMinutes','setUTCSeconds'];
function dateUI(){fields.forEach((f,i)=>{$(f).value=String(date[getters[i]]()+(i===1?1:0)).padStart(i===0?4:2,'0');});}
fields.forEach((f,i)=>{const wrap=document.createElement('div');wrap.className='datepart';wrap.innerHTML=`<input id="${f}" type="text" inputmode="numeric" aria-label="UTC ${labels[i]}"><div class="date-arrows"><button type="button" aria-label="${labels[i]} 증가">▲</button><button type="button" aria-label="${labels[i]} 감소">▼</button></div><span>${labels[i]}</span>`;$('datetime').append(wrap);
 function change(delta){date[setters[i]](date[getters[i]]()+delta);dateUI();schedule();}
 wrap.querySelectorAll('button').forEach((b,j)=>b.onclick=()=>change(j===0?1:-1));$(f).onkeydown=e=>{if(e.key==='ArrowUp'||e.key==='ArrowDown'){e.preventDefault();change(e.key==='ArrowUp'?1:-1);}};
 $(f).onchange=()=>{const v=Number($(f).value);if(Number.isInteger(v)){date[setters[i]](v-(i===1?1:0));}dateUI();schedule();};});
date.setUTCMilliseconds(0);dateUI();$('now').onclick=()=>{date=new Date();date.setUTCMilliseconds(0);dateUI();schedule();};
let timer;function schedule(){clearTimeout(timer);timer=setTimeout(calculate,250);}
function panels(){$('zodi-controls').hidden=$('kind').value==='airglow';$('dgl-hint').hidden=!['dgl','zodiacal_dgl'].includes($('kind').value);$('air-controls').hidden=$('kind').value!=='airglow';}
for(const id of ['kind','band','nside','anti','elongation','height','sza','model'])$(id).onchange=()=>{panels();schedule();};
$('contour').onchange=()=>{if(result)draw();};$('calculate').onclick=calculate;
function stop(){generation++;worker?.terminate();worker=null;$('cancel').disabled=true;$('progress').hidden=true;}
$('cancel').onclick=()=>{stop();$('status').textContent='계산을 중지했습니다. 표시 중인 지도는 이전 결과입니다.';};
function request(){const r={kind:$('kind').value,band:$('band').value,nside:Number($('nside').value),utc:date.toISOString(),anti:$('anti').checked,elongation:Number($('elongation').value),height:Number($('height').value),sza:Number($('sza').value),model:$('model').value};
 if(r.elongation<0||r.elongation>180||!Number.isFinite(r.elongation))throw Error('Solar elongation은 0–180°입니다.');if(r.kind==='airglow'&&(!(r.height>0&&r.height<1500)||!(r.sza>=0&&r.sza<=180)))throw Error('H는 0–1500 km 사이, SZA는 0–180°입니다.');return r;}
function calculate(){if(!manifest)return;try{const r=request();stop();const id=generation;worker=new Worker(new URL('./worker.mjs?v=dgl-selfabs-1',import.meta.url),{type:'module'});$('cancel').disabled=false;$('progress').hidden=false;$('progress').value=0;$('status').textContent='데이터 읽기 및 계산 중…';const start=performance.now();
 worker.onmessage=({data:d})=>{if(id!==generation)return;if(d.error){stop();$('status').textContent=d.error;return;}if(d.progress!==undefined){$('progress').value=d.progress;$('status').textContent=`계산 중 · ${Math.round(d.progress*100)}%`;return;}result=d.result;savedRequest=r;stop();draw();$('status').textContent=`완료 · ${result.flux.length.toLocaleString()} pixels · ${((performance.now()-start)/1000).toFixed(1)} s`;for(const k of ['fits','pdf','png'])$(k).disabled=false;};
 worker.onerror=e=>{stop();$('status').textContent='계산 실패: '+e.message;};worker.postMessage({request:r,manifest});}catch(e){$('status').textContent=e.message;}}
const area={x:90,y:52,w:1020,h:510},gridW=340,gridH=170;
function sample(v){const i=vectorPixel(savedRequest.nside,v);return {mag:result.mag[i],i};}
const stops=[[253,231,37],[94,201,98],[33,145,140],[59,82,139],[68,1,84]];
function color(v,lo,hi){const u=Math.max(0,Math.min(4,4*(v-lo)/(hi-lo))),i=Math.min(3,Math.floor(u)),f=u-i;return stops[i].map((x,j)=>Math.round(x*(1-f)+stops[i+1][j]*f));}
let scale;
function draw(){const r=savedRequest,b=manifest.bands[r.band];$('map-title').textContent=`${({zodiacal:'Zodiacal light',dgl:'DGL',zodiacal_dgl:'Zodiacal light + DGL',airglow:'Airglow'})[r.kind]} | ${r.band} · ${b[0]}–${b[1]} nm | ${r.kind!=='airglow'?'UTC '+r.utc.replace('T',' ').slice(0,19):'H = '+r.height+' km · SZA = '+r.sza+'° · mean Sun'} | NSIDE ${r.nside}${r.kind!=='airglow'&&r.anti?' · anti-sun ≥ '+r.elongation+'°':''}`;
 const finite=Array.from(result.mag).filter(Number.isFinite);let lo=Infinity,hi=-Infinity;for(const x of finite){lo=Math.min(lo,x);hi=Math.max(hi,x);}if(!finite.length){lo=18;hi=24;}lo=Math.floor(lo*2)/2;hi=Math.max(lo+.5,Math.ceil(hi*2)/2);scale={lo,hi};
 ctx.fillStyle='white';ctx.fillRect(0,0,canvas.width,canvas.height);const image=ctx.createImageData(area.w,area.h);
 for(let y=0;y<area.h;y++)for(let x=0;x<area.w;x++){const v=inverseMoll(2*(x+.5)/area.w-1,1-2*(y+.5)/area.h);if(!v)continue;const {mag,i}=sample(v);let rgb=Number.isFinite(mag)?color(mag,lo,hi):[245,245,245];if(result.earth?.[i])rgb=[174,181,190];else if(result.valid&&!result.valid[i]&&(x+y)%12<2)rgb=[110,110,110];const p=(y*area.w+x)*4;image.data.set([...rgb,255],p);}ctx.putImageData(image,area.x,area.y);
 ctx.strokeStyle='#ffffff60';ctx.lineWidth=1;
 function path(points){ctx.beginPath();let first=true;for(const [x,y]of points){const px=area.x+(x+1)*area.w/2,py=area.y+(1-y)*area.h/2;if(first){ctx.moveTo(px,py);first=false;}else ctx.lineTo(px,py);}ctx.stroke();}
 for(const dec of [-60,-30,0,30,60])path(Array.from({length:181},(_,i)=>moll((-180+i*2)*D,dec*D)));
 for(const ra of [-150,-120,-90,-60,-30,0,30,60,90,120,150])path(Array.from({length:181},(_,i)=>moll(ra*D,(-90+i)*D)));
 if($('contour').checked)contours(lo,hi);
 ctx.lineWidth=1;ctx.strokeStyle='#5b6b7a';ctx.beginPath();ctx.ellipse(area.x+area.w/2,area.y+area.h/2,area.w/2,area.h/2,0,0,2*Math.PI);ctx.stroke();
 ctx.fillStyle='#203348';ctx.font='14px system-ui';ctx.textAlign='center';for(const ra of [150,90,30,0,330,270,210]){const [x]=moll(ra*D,0);ctx.fillText(ra+'°',area.x+(x+1)*area.w/2,area.y+area.h+24);}ctx.fillText('RA (Equatorial) — increases to the left',600,area.y+area.h+48);ctx.save();ctx.translate(25,area.y+area.h/2);ctx.rotate(-Math.PI/2);ctx.fillText('DEC (°)',0,0);ctx.restore();for(const dec of [-60,-30,0,30,60]){const [,y]=moll(0,dec*D);ctx.textAlign='right';ctx.fillText(dec+'°',area.x-12,area.y+(1-y)*area.h/2+5);}
 for(let x=0;x<700;x++){ctx.fillStyle=`rgb(${color(lo+(hi-lo)*x/699,lo,hi)})`;ctx.fillRect(250+x,644,1,18);}ctx.fillStyle='#203348';ctx.textAlign='center';for(let i=0;i<=6;i++){const x=250+700*i/6;ctx.fillText((lo+(hi-lo)*i/6).toFixed(1),x,681);}ctx.fillText('AB mag / arcsec²',600,707);
 $('legend').textContent=`Contour: 0.5 mag · ${lo.toFixed(1)}–${hi.toFixed(1)} · 큰 AB mag = 어두운 색${r.kind==='airglow'?' · 회색: 지구 차폐 · 빗금: 모델 범위 밖 / photoelectron cutoff · 기준 Sun RA = 0°, DEC = 0°':''}`;}
function contours(lo,hi){const w=gridW,h=gridH,a=new Float32Array((w+1)*(h+1));a.fill(NaN);for(let y=0;y<=h;y++)for(let x=0;x<=w;x++){const v=inverseMoll(2*x/w-1,1-2*y/h);if(v){const s=sample(v);if(!result.valid||result.valid[s.i])a[y*(w+1)+x]=s.mag;}}
 const curves=[];
 for(let level=lo+.5;level<hi;level+=.5)for(const path of contourPaths(a,w,h,level,area))curves.push({...path,level});
 ctx.lineWidth=.8;ctx.strokeStyle='#ffffffb0';ctx.beginPath();
 for(const path of curves){ctx.moveTo(...path.points[0]);for(let i=1;i<path.points.length;i++)ctx.lineTo(...path.points[i]);}ctx.stroke();
 const dense=curves.length>200,labels=contourLabels(curves,area,dense?21:28,dense?11:14);
 // Draw callouts first, then all label boxes so leaders cannot cover the numbers.
 for(const label of labels)if(label.callout){
  const x=Math.max(label.rect.x,Math.min(label.rect.x+label.rect.w,label.anchor[0])),y=Math.max(label.rect.y,Math.min(label.rect.y+label.rect.h,label.anchor[1]));
  ctx.beginPath();ctx.moveTo(...label.anchor);ctx.lineTo(x,y);ctx.lineWidth=dense?1.4:2.5;ctx.strokeStyle='#203348';ctx.stroke();ctx.lineWidth=dense?.5:1;ctx.strokeStyle='white';ctx.stroke();
  ctx.beginPath();ctx.arc(...label.anchor,dense?.8:1.8,0,2*Math.PI);ctx.fillStyle='white';ctx.fill();
 }
 ctx.font=(dense?'8':'10')+'px system-ui';ctx.textAlign='center';ctx.textBaseline='middle';
 for(const label of labels){if(dense){ctx.strokeStyle='#203348';ctx.lineJoin='round';ctx.lineWidth=2.5;ctx.strokeText(label.level.toFixed(1),label.x,label.y);}else{ctx.fillStyle='#203348f0';ctx.fillRect(label.rect.x,label.rect.y,label.rect.w,label.rect.h);ctx.strokeStyle='#ffffffb0';ctx.lineWidth=.5;ctx.strokeRect(label.rect.x,label.rect.y,label.rect.w,label.rect.h);}ctx.fillStyle='white';ctx.fillText(label.level.toFixed(1),label.x,label.y);}
 ctx.textBaseline='alphabetic';
}

canvas.onmousemove=e=>{if(!result)return;const box=canvas.getBoundingClientRect(),x=(e.clientX-box.left)*canvas.width/box.width,y=(e.clientY-box.top)*canvas.height/box.height,v=inverseMoll(2*(x-area.x)/area.w-1,1-2*(y-area.y)/area.h);if(!v){$('cursor').textContent='mag at (RA, DEC) = (—, —)';return;}const {mag,i}=sample(v),ra=mod(Math.atan2(v[1],v[0])/D,360),dec=Math.asin(v[2])/D;$('cursor').textContent=`mag at (RA, DEC) = (${ra.toFixed(2)}°, ${dec.toFixed(2)}°) : ${Number.isFinite(mag)?mag.toFixed(3):mag===Infinity?'∞':'—'}${result.earth?.[i]?' · Earth occulted':result.valid&&!result.valid[i]?' · outside model support':''}`;};
function filename(){const r=savedRequest,stamp=r.kind!=='airglow'?r.utc.slice(0,19).replaceAll('-','').replaceAll(':','').replace('T','_'):`H${r.height}_SZA${r.sza}_${r.model}`;return `${r.kind}_${r.band}_nside${r.nside}_${stamp}${r.kind!=='airglow'&&r.anti?'_anti-sun':''}`;}
function download(blob,name){const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),10000);}
$('fits').onclick=async()=>{try{const bytes=fitsMap(result,savedRequest,manifest),blob=await new Response(new Blob([bytes]).stream().pipeThrough(new CompressionStream('gzip'))).blob();download(blob,filename()+'.fits.gz');}catch(e){$('status').textContent='FITS 저장 실패: '+e.message;}};
$('pdf').onclick=()=>{document.title=filename();window.print();};$('png').onclick=()=>canvas.toBlob(blob=>download(blob,filename()+'.png'));
try{if(!('Worker'in window)||!('DecompressionStream'in window)||!('CompressionStream'in window))throw Error('최신 Chrome, Edge, Firefox 또는 Safari가 필요합니다.');const response=await fetch(new URL('../data/manifest.json?v=dgl-selfabs-1',import.meta.url));if(!response.ok)throw Error('manifest 읽기 실패');manifest=await response.json();$('anti').checked=true;$('elongation').value='60';panels();calculate();}catch(e){$('status').textContent=e.message+' · HTML을 직접 열지 말고 HTTP 서버 또는 GitHub Pages에서 실행하십시오.';}
