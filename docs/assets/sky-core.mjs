// S-DRIFT browser routines, modified 2026-10-08. GLOW-derived integration:
// This software is part of the GLOW model. Use is governed by the Open Source
// Academic Research License Agreement contained in data/glowlicense.txt.
export const R=6371, D=Math.PI/180, ARC=206264.80624709636;
const clamp=(v,a,b)=>Math.max(a,Math.min(b,v));
export const mod=(x,n)=>((x%n)+n)%n;
export function ring(n,r){
  const np=12*n*n;
  if(r<n) return {start:2*r*(r-1),count:4*r,z:1-r*r/(3*n*n),offset:.5};
  if(r<=3*n) return {start:2*n*(n-1)+(r-n)*4*n,count:4*n,z:(2*n-r)*2/(3*n),offset:(r+n)%2?.0:.5};
  const k=4*n-r;return {start:np-2*k*(k+1),count:4*k,z:-1+k*k/(3*n*n),offset:.5};
}
export function pixVector(n,p){
  const cap=2*n*(n-1),np=12*n*n;let r;
  if(p<cap) r=Math.floor((1+Math.sqrt(1+2*p))/2);
  else if(p<np-cap) r=Math.floor((p-cap)/(4*n))+n;
  else r=4*n-Math.floor((1+Math.sqrt(2*(np-p)-1))/2);
  const row=ring(n,r),phi=(p-row.start+row.offset)*2*Math.PI/row.count,s=Math.sqrt(1-row.z*row.z);
  return [s*Math.cos(phi),s*Math.sin(phi),row.z];
}
export function vectorPixel(n,v){
  const z=v[2],za=Math.abs(z),tt=mod(Math.atan2(v[1],v[0]),2*Math.PI)/(Math.PI/2),cap=2*n*(n-1);
  if(za<=2/3){
    const jp=Math.floor(n*(.5+tt-.75*z)),jm=Math.floor(n*(.5+tt+.75*z));
    const ir=n+1+jp-jm,kshift=1-(ir&1),ip=mod(Math.floor((jp+jm-n+kshift+1)/2),4*n);
    return cap+(ir-1)*4*n+ip;
  }
  const tp=tt-Math.floor(tt),tmp=n*Math.sqrt(3*(1-za)),jp=Math.floor(tp*tmp),jm=Math.floor((1-tp)*tmp),ir=jp+jm+1;
  const ip=Math.min(4*ir-1,Math.floor(tt*ir));
  return z>0?2*ir*(ir-1)+ip:12*n*n-2*ir*(ir+1)+ip;
}
export function sampleMap(a,n,v){
  const theta=Math.acos(clamp(v[2],-1,1)),phi=mod(Math.atan2(v[1],v[0]),2*Math.PI);
  const sample=r=>{
    const row=ring(n,r),u=phi*row.count/(2*Math.PI)-row.offset,j=Math.floor(u),f=u-j;
    return a[row.start+mod(j,row.count)]*(1-f)+a[row.start+mod(j+1,row.count)]*f;
  };
  const first=Math.acos(ring(n,1).z),last=Math.PI-first;
  if(theta<first||theta>last){
    const r=theta<first?1:4*n-1,row=ring(n,r),f=(theta<first?theta:Math.PI-theta)/first;
    let mean=0;for(let i=0;i<4;i++)mean+=a[row.start+i]/4;
    return mean*(1-f)+sample(r)*f;
  }
  let lo=1,hi=4*n-1;
  while(hi-lo>1){const m=(lo+hi)>>1;if(Math.acos(ring(n,m).z)<=theta)lo=m;else hi=m;}
  const ta=Math.acos(ring(n,lo).z),tb=Math.acos(ring(n,hi).z),f=(theta-ta)/(tb-ta);
  return sample(lo)*(1-f)+sample(hi)*f;
}
export function rotation(v,angle){
  const e=23.439291111*D,c=Math.cos(e),s=Math.sin(e),y=c*v[1]+s*v[2],z=-s*v[1]+c*v[2];
  const x2=Math.cos(angle)*v[0]-Math.sin(angle)*y,y2=Math.sin(angle)*v[0]+Math.cos(angle)*y;
  return [x2,c*y2-s*z,s*y2+c*z];
}
export function sunLongitude(v){const e=23.439291111*D;return Math.atan2(Math.cos(e)*v[1]+Math.sin(e)*v[2],v[0]);}
export function solarAt(unix,data,meta){
  const u=(unix-meta.start_unix)/meta.step_seconds,i=Math.floor(u),f=u-i;
  if(i<0||i>=meta.count-1)throw Error('날짜는 1900–2100년 범위에서 선택하십시오.');
  const v=[0,1,2].map(j=>data[i*3+j]*(1-f)+data[(i+1)*3+j]*f),norm=Math.hypot(...v);
  return v.map(x=>x/norm);
}
export function abZodi(value){return value>0?-2.5*Math.log10(value*1e6/ARC**2/3631):Infinity;}
export function abAir(value,band){return value>0?-2.5*Math.log10(6.62607015e-27*value/(Math.log(band[1]/band[0])*ARC**2)/(3631e-23)):Infinity;}
export function geometry(H,sza,v,top){
  const c=Math.cos(sza*D),s=Math.sin(sza*D),rr=R+H,b=rr*(c*v[0]+s*v[2]);
  const occulted=b<0&&b*b-rr*rr+R*R>=0,hmin=occulted?0:Math.sqrt(Math.max(rr*rr-Math.min(b,0)**2,0))-R;
  const disc=b*b-rr*rr+(R+top)**2,root=Math.sqrt(Math.max(0,disc));
  return {occulted,hmin,start:Math.max(-b-root,0),end:disc>=0?Math.max(-b+root,0):0,p:[rr*c,0,rr*s]};
}
export function bracket(axis,x){let lo=0,hi=axis.length-1;while(hi-lo>1){const m=(lo+hi)>>1;if(axis[m]<=x)lo=m;else hi=m;}return lo;}
export function integrateRay(H,sza,v,q,info,step=2){
  const z=info.altitude_km,chi=info.chi_deg,nz=z.length,g=geometry(H,sza,v,z.at(-1));
  if(g.occulted)return {flux:NaN,lit:NaN,valid:false,cutoff:false,hmin:0,earth:true};
  const count=Math.max(1,Math.ceil((g.end-g.start)/step)),ds=(g.end-g.start)/count,factor=ds*1e5/(4*Math.PI);
  let flux=0,lit=0,cutoff=false;
  for(let j=0;j<count;j++){
    const t=g.start+ds*(j+.5),x=g.p[0]+t*v[0],y=t*v[1],zz=g.p[2]+t*v[2],rr=Math.hypot(x,y,zz),h=rr-R;
    if(h<z[0]||h>z.at(-1))continue;
    const angle=Math.acos(clamp(x/rr,-1,1))/D,shadow=x<0&&y*y+zz*zz<=R*R;
    if(!shadow&&angle>=1.85/D)cutoff=true;
    const iz=bracket(z,h),ic=bracket(chi,angle),fz=(h-z[iz])/(z[iz+1]-z[iz]),fc=(angle-chi[ic])/(chi[ic+1]-chi[ic]);
    const a=(ic*nz+iz)*2,b=a+2,c=a+nz*2,d=c+2;
    let value=(1-fc)*((1-fz)*q[a]+fz*q[b])+fc*((1-fz)*q[c]+fz*q[d]);
    if(shadow)value+=(1-fc)*((1-fz)*q[a+1]+fz*q[b+1])+fc*((1-fz)*q[c+1]+fz*q[d+1]);
    flux+=value;if(!shadow)lit+=value;
  }
  return {flux:flux*factor,lit:lit*factor,earth:false,hmin:g.hmin,cutoff,
    valid:g.hmin>=(info.metadata.support_min_altitude_km??z[0])&&!cutoff};
}
export function inverseMoll(x,y){
  if(x*x+y*y>1)return null;
  const t=Math.asin(y),dec=Math.asin(clamp((2*t+Math.sin(2*t))/Math.PI,-1,1)),lon=Math.PI*x/Math.max(1e-12,Math.cos(t));
  return [Math.cos(dec)*Math.cos(-lon),Math.cos(dec)*Math.sin(-lon),Math.sin(dec)];
}
export function moll(ra,dec){
  let t=dec;for(let i=0;i<12;i++){const den=2+2*Math.cos(2*t);if(den<1e-10)break;t-=(2*t+Math.sin(2*t)-Math.PI*Math.sin(dec))/den;}
  return [-mod(ra+Math.PI,2*Math.PI)+Math.PI, t].map((a,i)=>i?Math.sin(a):a*Math.cos(t)/Math.PI);
}

// FITS uses big-endian IEEE floats and 2880-byte blocks. All sky exports are
// single-pixel-per-row RING HEALPix BINTABLEs, independently checked by Astropy.
function card(key,val){
  let text=key.padEnd(8)+'= '+(typeof val==='string'?"'"+val.replaceAll("'","''")+"'":typeof val==='boolean'?(val?'T':'F').padStart(20):String(val).padStart(20));
  return text.padEnd(80).slice(0,80);
}
function header(rows){const text=[...rows.map(([k,v])=>card(k,v)),'END'.padEnd(80)].join('');return new TextEncoder().encode(text.padEnd(Math.ceil(text.length/2880)*2880));}
export function fitsMap(result,request,manifest){
  const n=result.flux.length,air=request.kind==='airglow',hasDgl=['dgl','zodiacal_dgl'].includes(request.kind),hasZodi=['zodiacal','zodiacal_dgl'].includes(request.kind);
  const names=air?['AB_MAG_ARCSEC2','PHOTON_RADIANCE','SUNLIT_RADIANCE','EARTH_OCCULTED','MODEL_VALID','MIN_ALTITUDE','PHOTO_CUTOFF']:['AB_MAG_ARCSEC2','INTENSITY_MJY_SR'];
  if(request.kind==='zodiacal_dgl')names.push('ZODIACAL_MJY_SR','DGL_MJY_SR');
  const primary=header([['SIMPLE',true],['BITPIX',8],['NAXIS',0],['EXTEND',true]]);
  const rows=[['XTENSION','BINTABLE'],['BITPIX',8],['NAXIS',2],['NAXIS1',names.length*4],['NAXIS2',n],['PCOUNT',0],['GCOUNT',1],['TFIELDS',names.length],['PIXTYPE','HEALPIX'],['ORDERING','RING'],['NSIDE',request.nside],['FIRSTPIX',0],['LASTPIX',n-1],['INDXSCHM','IMPLICIT'],['COORDSYS','C'],['BAND',request.band],['WEBVER',2],['CREATOR','S-DRIFT HTML'],['DATATYPE',request.kind.toUpperCase()],['BANDMIN',manifest.bands[request.band][0]/1000],['BANDMAX',manifest.bands[request.band][1]/1000],['PHOTCNT',true]];
  names.forEach((name,i)=>{rows.push(['TTYPE'+(i+1),name],['TFORM'+(i+1),'1E']);});
  rows.push(['TUNIT1','mag/arcsec2'],['TUNIT2',air?'photons/s/cm2/sr':'MJy/sr']);
  if(air){const m=manifest.airglow[request.model];rows.push(['SAT_H',request.height],['SAT_SZA',request.sza],['DS_KM',2],['INTVER',2],['FULLBAND',false],['ABSORPT',false],['TWILHIST',false],['NIGHTAPP',request.model==='combined'],['TABLESHA',m.source_sha256],['MODEL',m.metadata.model],['TUNIT3','photons/s/cm2/sr'],['TUNIT6','km']);}
  else{
    rows.push(['DATE-UTC',request.utc.replace('Z','')],['ANTISUN',request.anti],['ELONGMIN',request.elongation]);
    if(hasZodi)rows.push(['REF-UTC',manifest.solar.reference_utc],['ROT-DEG',result.delta/D],['SRC_SHA',manifest.zodiacal[request.band].source_sha256]);
    if(hasDgl){const m=manifest.dgl;rows.push(['CORRMOD','C2022_all'],['CORRSCL',m.metadata.scale],['CORRMIN',m.metadata.correlation_min_nm],['UVEXT','CONSTANT'],['BETAOP','DIVIDE'],['DGLFIXED',true],['DGL_SHA',m.bands[request.band].source_sha256]);}
    if(request.kind==='zodiacal_dgl')rows.push(['TUNIT3','MJy/sr'],['TUNIT4','MJy/sr']);
  }
  const ext=header(rows),dataBytes=n*names.length*4,out=new Uint8Array(primary.length+ext.length+Math.ceil(dataBytes/2880)*2880);
  out.set(primary);out.set(ext,primary.length);const view=new DataView(out.buffer,primary.length+ext.length);
  for(let i=0;i<n;i++){
    const vals=air?[result.mag[i],result.flux[i],result.lit[i],result.earth[i],result.valid[i],result.hmin[i],result.cutoff[i]]:[result.mag[i],result.flux[i]];
    if(request.kind==='zodiacal_dgl')vals.push(result.zodiacal[i],result.dgl[i]);
    for(let j=0;j<vals.length;j++)view.setFloat32((i*vals.length+j)*4,vals[j],false);
  }
  return out;
}
