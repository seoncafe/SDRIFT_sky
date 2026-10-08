// Modified 2026-10-08. See sky-core.mjs and ../data/glowlicense.txt.
import {pixVector,rotation,solarAt,sunLongitude,abZodi,abAir,integrateRay,D,mod} from './sky-core.mjs';
const cache=new Map();
async function array(meta){
 if(!cache.has(meta.file))cache.set(meta.file,(async()=>{
 const r=await fetch(new URL('../data/'+meta.file,import.meta.url));if(!r.ok)throw Error(`데이터 읽기 실패: ${r.status} ${meta.file}`);
 const b=await new Response(r.body.pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
 if(b.byteLength!==meta.bytes)throw Error('데이터 크기가 일치하지 않습니다.');
 const hash=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',b)),x=>x.toString(16).padStart(2,'0')).join('');
 if(hash!==meta.sha256)throw Error('데이터 SHA-256 검증 실패');
 return new Float32Array(b);})());return cache.get(meta.file);
}
self.onmessage=async({data:{request:r,manifest:m}})=>{try{
 const np=12*r.nside*r.nside,result={flux:new Float32Array(np),mag:new Float32Array(np)};
 let q,info,sun,delta;
 if(r.kind==='zodiacal'){q=await array(m.zodiacal[r.band]);sun=solarAt(Date.parse(r.utc)/1000,await array(m.solar),m.solar);delta=mod(sunLongitude(sun)-sunLongitude(m.solar.reference_vector)+Math.PI,2*Math.PI)-Math.PI;result.delta=delta;result.sun=sun;}
 else {info=m.airglow[r.model];q=await array(info.bands[r.band]);for(const key of ['lit','hmin'])result[key]=new Float32Array(np);for(const key of ['earth','valid','cutoff'])result[key]=new Uint8Array(np);}
 const {sampleMap}=await import('./sky-core.mjs');
 for(let i=0;i<np;i++){
 const v=pixVector(r.nside,i);
 if(r.kind==='zodiacal'){
 const masked=r.anti&&v.reduce((s,x,j)=>s+x*sun[j],0)>Math.cos(r.elongation*D);
 result.flux[i]=masked?NaN:sampleMap(q,128,rotation(v,-delta));result.mag[i]=abZodi(result.flux[i]);if(masked)result.mag[i]=NaN;
 }else{const a=integrateRay(r.height,r.sza,v,q,info);for(const key of ['flux','lit','hmin','earth','valid','cutoff'])result[key][i]=a[key];result.mag[i]=a.earth?NaN:abAir(a.flux,m.bands[r.band]);}
 if(i%1024===0){self.postMessage({progress:i/np});await new Promise(resolve=>setTimeout(resolve,0));}
 }
 self.postMessage({result},Object.values(result).filter(x=>ArrayBuffer.isView(x)).map(x=>x.buffer));
 }catch(e){self.postMessage({error:e.message});}};
