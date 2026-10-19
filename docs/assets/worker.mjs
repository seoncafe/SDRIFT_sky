// Modified 2026-10-20. See sky-core.mjs and ../data/glowlicense.txt.
import {pixVector,rotation,solarAt,bodyVector,topocentricMoon,sunLongitude,zodiacalTemplate,abZodi,abAir,integrateRay,orbitState,referenceDirection,visibility,geometry,airglowMjy,toEquatorial,D,mod} from './sky-core.mjs?v=sso-harmonic-3';
const cache=new Map();
async function array(meta){
 if(!cache.has(meta.file))cache.set(meta.file,(async()=>{
 const r=await fetch(new URL('../data/'+meta.file+'?sha='+meta.sha256,import.meta.url));if(!r.ok)throw Error(`데이터 읽기 실패: ${r.status} ${meta.file}`);
 const b=await new Response(r.body.pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
 if(b.byteLength!==meta.bytes)throw Error('데이터 크기가 일치하지 않습니다.');
 const hash=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',b)),x=>x.toString(16).padStart(2,'0')).join('');
 if(hash!==meta.sha256)throw Error('데이터 SHA-256 검증 실패');
 return new Float32Array(b);})().catch(e=>{cache.delete(meta.file);throw e;}));return cache.get(meta.file);
}
// Requests carry an id; a newer request or {cancel:true} makes the running one stop at its next yield.
let latest;
self.onmessage=async({data})=>{if(data.cancel){latest=Symbol('cancelled');return;}
 const {id,request:r,manifest:m}=data;latest=id;try{
 const np=12*r.nside*r.nside,result={flux:new Float32Array(np),mag:new Float32Array(np)};
 const air=['airglow','total'].includes(r.kind),total=r.kind==='total';
 let q,info,sun,delta,dgl,moon,state,zodi;
 sun=solarAt(Date.parse(r.utc)/1000,await array(m.solar),m.solar);result.sun=sun;
 state=orbitState(sun,r.height,r.orbit,r.orbitalPhaseDeg,Date.parse(r.utc)/1000,r.evolve?Date.parse(r.epoch)/1000:null);result.state=state;
 moon=topocentricMoon(bodyVector(Date.parse(r.utc)/1000,await array(m.lunar),m.lunar),state.zenith,r.height);result.moon=moon;
 for(const key of ['hmin','theta','elv'])result[key]=new Float32Array(np);
 for(const key of ['earth','valid','view'])result[key]=new Uint8Array(np);
 if(!air){
   if(r.kind==='zodiacal'||r.kind==='zodiacal_dgl')q=zodiacalTemplate(await array(m.zodiacal[r.band]),12*m.zodiacal[r.band].nside**2,sunLongitude(sun));
   if(r.kind==='dgl'||r.kind==='zodiacal_dgl')dgl=await array(m.dgl.bands[r.band]);
   delta=q?mod(sunLongitude(sun)-sunLongitude(m.solar.reference_vector)+Math.PI,2*Math.PI)-Math.PI:0;
   result.delta=delta;
   if(r.kind==='zodiacal_dgl'){result.zodiacal=new Float32Array(np);result.dgl=new Float32Array(np);}
 }else {info=m.airglow[r.model];q=await array(info.bands[r.band]);result.lit=new Float32Array(np);result.cutoff=new Uint8Array(np);}
 if(total){
  zodi=zodiacalTemplate(await array(m.zodiacal[r.band]),12*m.zodiacal[r.band].nside**2,sunLongitude(sun));dgl=await array(m.dgl.bands[r.band]);
  delta=mod(sunLongitude(sun)-sunLongitude(m.solar.reference_vector)+Math.PI,2*Math.PI)-Math.PI;result.delta=delta;
  for(const key of ['zodiacal','dgl','airglow','airglowPhoton']){result[key]=new Float32Array(np);result[key].fill(NaN);}
 }
 const {sampleMap}=await import('./sky-core.mjs?v=sso-harmonic-3');
 for(let i=0;i<np;i++){
 const v=toEquatorial(pixVector(r.nside,i),r.coords);
 const mask=visibility(v,state,r,moon);result.theta[i]=mask.theta;result.elv[i]=mask.elv;result.view[i]=mask.view;
 if(!air){
 const g0=geometry(r.height,state.sza,referenceDirection(v,state),1500);
 result.earth[i]=g0.occulted;result.valid[i]=!g0.occulted;result.hmin[i]=g0.hmin;
 const masked=!mask.view||g0.occulted;
 const z=q?sampleMap(q,m.zodiacal[r.band].nside,rotation(v,-delta)):0,g=dgl?sampleMap(dgl,m.dgl.bands[r.band].nside,v):0;
 result.flux[i]=masked?NaN:z+g;result.mag[i]=masked?NaN:abZodi(result.flux[i]);
 if(r.kind==='zodiacal_dgl'){result.zodiacal[i]=masked?NaN:z;result.dgl[i]=masked?NaN:g;}
 }else{
 result.flux[i]=NaN;result.mag[i]=NaN;result.lit[i]=NaN;
 if(!mask.view){const g0=geometry(r.height,state.sza,referenceDirection(v,state),1500);result.earth[i]=g0.occulted;result.hmin[i]=g0.hmin;}
 else{const a=integrateRay(r.height,state.sza,referenceDirection(v,state),q,info);for(const key of ['flux','lit','hmin','earth','valid','cutoff'])result[key][i]=a[key];result.mag[i]=a.earth?NaN:abAir(a.flux,m.bands[r.band]);
 if(total&&!a.earth){
  const z=sampleMap(zodi,m.zodiacal[r.band].nside,rotation(v,-delta)),g=sampleMap(dgl,m.dgl.bands[r.band].nside,v),ag=airglowMjy(a.flux,m.bands[r.band]);
  result.airglowPhoton[i]=a.flux;result.airglow[i]=ag;result.zodiacal[i]=z;result.dgl[i]=g;
  result.flux[i]=z+g+ag;result.mag[i]=abZodi(z+g+ag);
 }}
 }
 if(i%1024===0){self.postMessage({id,progress:i/np});await new Promise(resolve=>setTimeout(resolve,0));if(latest!==id)return;}
 }
 if(latest!==id)return;
 self.postMessage({id,result},Object.values(result).filter(x=>ArrayBuffer.isView(x)).map(x=>x.buffer));
 }catch(e){if(latest===id)self.postMessage({id,error:e.message});}};
