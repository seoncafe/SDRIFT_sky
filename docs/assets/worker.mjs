// Modified 2026-10-17. See sky-core.mjs and ../data/glowlicense.txt.
import {pixVector,rotation,solarAt,sunLongitude,abZodi,abAir,integrateRay,orbitState,referenceDirection,visibility,D,mod} from './sky-core.mjs?v=north-phase-1';
const cache=new Map();
async function array(meta){
 if(!cache.has(meta.file))cache.set(meta.file,(async()=>{
 const r=await fetch(new URL('../data/'+meta.file+'?sha='+meta.sha256,import.meta.url));if(!r.ok)throw Error(`데이터 읽기 실패: ${r.status} ${meta.file}`);
 const b=await new Response(r.body.pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
 if(b.byteLength!==meta.bytes)throw Error('데이터 크기가 일치하지 않습니다.');
 const hash=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',b)),x=>x.toString(16).padStart(2,'0')).join('');
 if(hash!==meta.sha256)throw Error('데이터 SHA-256 검증 실패');
 return new Float32Array(b);})());return cache.get(meta.file);
}
self.onmessage=async({data:{request:r,manifest:m}})=>{try{
 const np=12*r.nside*r.nside,result={flux:new Float32Array(np),mag:new Float32Array(np)};
 let q,info,sun,delta,dgl,moon,state;
 sun=solarAt(Date.parse(r.utc)/1000,await array(m.solar),m.solar);result.sun=sun;
 moon=solarAt(Date.parse(r.utc)/1000,await array(m.lunar),m.lunar);result.moon=moon;
 if(r.kind!=='airglow'){
   if(r.kind==='zodiacal'||r.kind==='zodiacal_dgl')q=await array(m.zodiacal[r.band]);
   if(r.kind==='dgl'||r.kind==='zodiacal_dgl')dgl=await array(m.dgl.bands[r.band]);
   sun=solarAt(Date.parse(r.utc)/1000,await array(m.solar),m.solar);
   delta=q?mod(sunLongitude(sun)-sunLongitude(m.solar.reference_vector)+Math.PI,2*Math.PI)-Math.PI:0;
   result.delta=delta;result.sun=sun;
   moon=solarAt(Date.parse(r.utc)/1000,await array(m.lunar),m.lunar);result.moon=moon;
   if(r.kind==='zodiacal_dgl'){result.zodiacal=new Float32Array(np);result.dgl=new Float32Array(np);}
 }else {info=m.airglow[r.model];q=await array(info.bands[r.band]);state=orbitState(sun,r.height,r.orbit,r.orbitalPhaseDeg,Date.parse(r.utc)/1000,r.evolve?Date.parse(r.epoch)/1000:null);result.state=state;for(const key of ['lit','hmin','theta'])result[key]=new Float32Array(np);for(const key of ['earth','valid','cutoff','view'])result[key]=new Uint8Array(np);}
 const {sampleMap}=await import('./sky-core.mjs?v=north-phase-1');
 for(let i=0;i<np;i++){
 const v=pixVector(r.nside,i);
 if(r.kind!=='airglow'){
 const masked=(r.anti&&v.reduce((s,x,j)=>s+x*sun[j],0)>Math.cos(r.elongation*D))||(moon&&r.moonAvoid&&v.reduce((s,x,j)=>s+x*moon[j],0)>Math.cos(r.moonAngle*D));
 const z=q?sampleMap(q,128,rotation(v,-delta)):0,g=dgl?sampleMap(dgl,128,v):0;
 result.flux[i]=masked?NaN:z+g;result.mag[i]=masked?NaN:abZodi(result.flux[i]);
 if(r.kind==='zodiacal_dgl'){result.zodiacal[i]=masked?NaN:z;result.dgl[i]=masked?NaN:g;}
 }else{
 const mask=visibility(v,state,r,moon);result.theta[i]=mask.theta;result.view[i]=mask.view;
 result.flux[i]=NaN;result.mag[i]=NaN;result.lit[i]=NaN;result.hmin[i]=NaN;
 if(mask.view){const a=integrateRay(r.height,state.sza,referenceDirection(v,state),q,info);for(const key of ['flux','lit','hmin','earth','valid','cutoff'])result[key][i]=a[key];result.mag[i]=a.earth?NaN:abAir(a.flux,m.bands[r.band]);}
 }
 if(i%1024===0){self.postMessage({progress:i/np});await new Promise(resolve=>setTimeout(resolve,0));}
 }
 self.postMessage({result},Object.values(result).filter(x=>ArrayBuffer.isView(x)).map(x=>x.buffer));
 }catch(e){self.postMessage({error:e.message});}};
