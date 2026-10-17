// S-DRIFT browser routines, modified 2026-10-17. GLOW-derived integration:
// This software is based on calculations using the GLOW model, governed by the Open Source
// Academic Research License Agreement contained in data/glowlicense.txt.
export const R=6371, D=Math.PI/180, ARC=206264.80624709636;
const clamp=(v,a,b)=>Math.max(a,Math.min(b,v));
export const mod=(x,n)=>((x%n)+n)%n;
// Same Astropy ICRS-axis -> Galactic basis rotation as sky_coordinates.py.
const EQ_TO_GAL=[[-0.054875657712591654,-0.8734370519556157,-0.48383507361671546],[0.49410943719272676,-0.4448297212232951,0.7469821839866675],[-0.8676661375596576,-0.19807633727300059,0.4559838136873016]];
export function fromEquatorial(v,coords='equatorial'){return coords==='galactic'?EQ_TO_GAL.map(a=>a.reduce((s,x,j)=>s+x*v[j],0)):v;}
export function toEquatorial(v,coords='equatorial'){return coords==='galactic'?[0,1,2].map(j=>EQ_TO_GAL.reduce((s,a,i)=>s+a[j]*v[i],0)):v;}

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
export function airglowMjy(value,band){return 6.62607015e-10*value/Math.log(band[1]/band[0]);}
export function abAir(value,band){return value>0?-2.5*Math.log10(6.62607015e-27*value/(Math.log(band[1]/band[0])*ARC**2)/(3631e-23)):Infinity;}

// Ideal solar-following circles for initial design; not an SGP4/J2 propagator.
export function orbitState(sun,H,orbit='dawn-dusk',orbitalPhaseDeg=0,unix=0,epoch=null){
 if(!['dawn-dusk','noon-midnight'].includes(orbit))throw Error('Unknown orbit');
 if(!Number.isFinite(H)||H<=0||H>=1500||!Number.isFinite(orbitalPhaseDeg))throw Error('Invalid orbit input');
 const dot=(a,b)=>a.reduce((v,x,j)=>v+x*b[j],0),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],unit=a=>{const n=Math.hypot(...a);return a.map(x=>x/n);};
 const s=unit(sun),k=Math.abs(s[2])<.999999?[0,0,1]:[1,0,0];
 const north=unit(k.map((x,j)=>x-dot(k,s)*s[j])),east=cross(s,north);
 const period=2*Math.PI*Math.sqrt((R+H)**3/398600.4418);
 const phi=mod(orbitalPhaseDeg+(epoch===null?0:(unix-epoch)*360/period),360),a=phi*D;
 const zenith=unit(s.map((x,j)=>Math.cos(a)*north[j]+Math.sin(a)*(orbit==='dawn-dusk'?east[j]:-x)));
 let c=clamp(dot(zenith,s),-1,1);if(Math.abs(Math.abs(c)-1)<1e-14)c=Math.sign(c);const sz=Math.acos(c)/D;
 const proj=zenith.map((x,j)=>x-c*s[j]),z=Math.hypot(...proj)>1e-10?unit(proj):north;
 return {sun:s,zenith,sza:sz,orbitalPhaseDeg:phi,period,limbZenithDeg:180-Math.asin(R/(R+H))/D,axes:[s,cross(z,s),z]};
}
export function referenceDirection(v,state){return state.axes.map(a=>a.reduce((s,x,j)=>s+x*v[j],0));}
export function visibility(v,state,r,moon){
 const dot=(a,b)=>a.reduce((s,x,j)=>s+x*b[j],0);
 const theta=Math.acos(clamp(dot(v,state.zenith),-1,1))/D;
 const elv=state.limbZenithDeg-theta;
 return {theta,elv,view:elv>=(r.minElvDeg??20)-1e-9&&
 (!r.anti||dot(v,state.sun)<=Math.cos(r.elongation*D)+1e-12)&&
 (!r.moonAvoid||dot(v,moon)<=Math.cos(r.moonAngle*D)+1e-12)};
}

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
  let text=key.padEnd(8)+'= '+(typeof val==='string'?"'"+val.replaceAll("'","''")+"'":typeof val==='boolean'?(val?'T':'F').padStart(20):String(val).toUpperCase().padStart(20));
  return text.padEnd(80).slice(0,80);
}
function header(rows){const text=[...rows.map(([k,v])=>card(k,v)),'END'.padEnd(80)].join('');return new TextEncoder().encode(text.padEnd(Math.ceil(text.length/2880)*2880));}
export function fitsMap(result,request,manifest){
  const n=result.flux.length,air=['airglow','total'].includes(request.kind),total=request.kind==='total',hasDgl=['dgl','zodiacal_dgl'].includes(request.kind),hasZodi=['zodiacal','zodiacal_dgl'].includes(request.kind);
  const names=air?['AB_MAG_ARCSEC2','PHOTON_RADIANCE','SUNLIT_RADIANCE','EARTH_OCCULTED','MODEL_VALID','MIN_ALTITUDE','PHOTO_CUTOFF','ZENITH_ANGLE','VIEW_ALLOWED','ELV']:['AB_MAG_ARCSEC2','INTENSITY_MJY_SR'];
  if(request.kind==='zodiacal_dgl')names.push('ZODIACAL_MJY_SR','DGL_MJY_SR');
  if(!air&&result.state)names.push('EARTH_OCCULTED','MODEL_VALID','MIN_ALTITUDE','ZENITH_ANGLE','VIEW_ALLOWED','ELV');
  if(total){names[1]='INTENSITY_MJY_SR';names.push('ZODIACAL_MJY_SR','DGL_MJY_SR','AIRGLOW_MJY_SR','AIRGLOW_PHOTON_RADIANCE');}
  const primary=header([['SIMPLE',true],['BITPIX',8],['NAXIS',0],['EXTEND',true]]);
  const rows=[['XTENSION','BINTABLE'],['BITPIX',8],['NAXIS',2],['NAXIS1',names.length*4],['NAXIS2',n],['PCOUNT',0],['GCOUNT',1],['TFIELDS',names.length],['PIXTYPE','HEALPIX'],['ORDERING','RING'],['NSIDE',request.nside],['FIRSTPIX',0],['LASTPIX',n-1],['INDXSCHM','IMPLICIT'],['COORDSYS',request.coords==='galactic'?'G':'C'],['BAND',request.band],['WEBVER',2],['CREATOR','S-DRIFT HTML'],['DATATYPE',request.kind.toUpperCase()],['BANDMIN',manifest.bands[request.band][0]/1000],['BANDMAX',manifest.bands[request.band][1]/1000],['PHOTCNT',true]];
  names.forEach((name,i)=>{rows.push(['TTYPE'+(i+1),name],['TFORM'+(i+1),'1E']);});
  rows.push(['TUNIT1','mag/arcsec2'],['TUNIT2',air&&!total?'photons/s/cm2/sr':'MJy/sr']);
  if(air){const m=manifest.airglow[request.model];rows.push(['SAT_H',request.height],['SAT_SZA',result.state?.sza??request.sza],['DS_KM',2],['INTVER',2],['FULLBAND',false],['ABSORPT',false],['TWILHIST',false],['NIGHTAPP',request.model==='combined'],['TABLESHA',m.source_sha256],['MODEL',m.metadata.model],['TUNIT3','photons/s/cm2/sr'],['TUNIT6','km'],['TUNIT8','deg'],['TUNIT10','deg'],['ELVMIN',request.minElvDeg??20],['ELVREF','SPHERICAL EARTH'],['EARTHR',R],['DATE-UTC',request.utc.replace('Z','')],['ORBIT',request.orbit],['ORBMODEL','IDEAL SOLAR CIRCLE'],['PHASE',result.state.orbitalPhaseDeg],['PHREF','NORTH'],['PHCONV','NORTH_V1'],['PHASE0',request.orbitalPhaseDeg],['PHMODE',request.evolve?'EPOCH':'MANUAL'],['EPOCH',request.epoch],['PERIOD',result.state.period],['SUNMASK',request.anti],['ELONGMIN',request.elongation],['MOONMASK',request.moonAvoid],['MOONMIN',request.moonAngle],['MOONOBS','GEOCENTRIC']);
    for(const [v,prefix]of [[result.sun,'SUN'],[result.moon,'MOON'],[result.state.zenith,'ZEN']])if(v)rows.push([prefix+'RA',mod(Math.atan2(v[1],v[0])/D,360)],[prefix+'DEC',Math.asin(clamp(v[2],-1,1))/D]);
  }
  else{
    rows.push(['DATE-UTC',request.utc.replace('Z','')],['ANTISUN',request.anti],['ELONGMIN',request.elongation]);
    if(result.sun)rows.push(['SUNRA',mod(Math.atan2(result.sun[1],result.sun[0])/D,360)],['SUNDEC',Math.asin(result.sun[2])/D]);
    if(result.moon)rows.push(['MOONMASK',!!request.moonAvoid],['MOONMIN',request.moonAngle??60],['MOONRA',mod(Math.atan2(result.moon[1],result.moon[0])/D,360)],['MOONDEC',Math.asin(result.moon[2])/D],['MOONOBS','GEOCENTRIC']);
    if(hasZodi)rows.push(['REF-UTC',manifest.solar.reference_utc],['ROT-DEG',result.delta/D],['SRC_SHA',manifest.zodiacal[request.band].source_sha256]);
    if(hasDgl){const m=manifest.dgl;rows.push(['CORRMOD','C2022_all'],['CORRSCL',m.metadata.scale],['CORRMIN',m.metadata.correlation_min_nm],['UVEXT','CONSTANT'],['BETAOP','MULTIPLY'],['ALPHAPP',m.metadata.alpha_approximation],['BETAREF',m.metadata.beta_reference],['DGLFIXED',true],['DGL_SHA',m.bands[request.band].source_sha256]);}
    if(request.kind==='zodiacal_dgl')rows.push(['TUNIT3','MJy/sr'],['TUNIT4','MJy/sr']);
  }
  if(total){
    const m=manifest.dgl;
    rows.push(['TUNIT11','MJy/sr'],['TUNIT12','MJy/sr'],['TUNIT13','MJy/sr'],['TUNIT14','photons/s/cm2/sr'],
      ['REF-UTC',manifest.solar.reference_utc],['ROT-DEG',result.delta/D],['SRC_SHA',manifest.zodiacal[request.band].source_sha256],
      ['DGL_SHA',m.bands[request.band].source_sha256],['CORRMOD','C2022_all'],['CORRSCL',m.metadata.scale],['BETAOP','MULTIPLY'],
      ['ALPHAPP',m.metadata.alpha_approximation],['BETAREF',m.metadata.beta_reference],['DGLFIXED',true]);
  }
  if(!air&&result.state){
    const s=result.state,v=s.zenith;
    rows.push(['ORBIT',request.orbit],['ORBMODEL','IDEAL SOLAR CIRCLE'],['SAT_H',request.height],['SAT_SZA',s.sza],
      ['PHASE',s.orbitalPhaseDeg],['PHASE0',request.orbitalPhaseDeg],['PHREF','NORTH'],['PHCONV','NORTH_V1'],
      ['PHMODE',request.evolve?'EPOCH':'MANUAL'],['PERIOD',s.period],['ELVMIN',request.minElvDeg??20],['ELVREF','SPHERICAL EARTH'],['EARTHR',R],
      ['ZENRA',mod(Math.atan2(v[1],v[0])/D,360)],['ZENDEC',Math.asin(clamp(v[2],-1,1))/D]);
    if(request.evolve)rows.push(['EPOCH',request.epoch]);
    for(const [name,unit]of [['MIN_ALTITUDE','km'],['ZENITH_ANGLE','deg'],['ELV','deg']])rows.push(['TUNIT'+(names.indexOf(name)+1),unit]);
  }
  const ext=header(rows),dataBytes=n*names.length*4,out=new Uint8Array(primary.length+ext.length+Math.ceil(dataBytes/2880)*2880);
  out.set(primary);out.set(ext,primary.length);const view=new DataView(out.buffer,primary.length+ext.length);
  for(let i=0;i<n;i++){
    const vals=air?[result.mag[i],result.flux[i],result.lit[i],result.earth[i],result.valid[i],result.hmin[i],result.cutoff[i],result.theta[i],result.view[i],result.elv[i]]:[result.mag[i],result.flux[i]];
    if(request.kind==='zodiacal_dgl')vals.push(result.zodiacal[i],result.dgl[i]);
    if(!air&&result.state)vals.push(result.earth[i],result.valid[i],result.hmin[i],result.theta[i],result.view[i],result.elv[i]);
    if(total)vals.push(result.zodiacal[i],result.dgl[i],result.airglow[i],result.airglowPhoton[i]);
    for(let j=0;j<vals.length;j++)view.setFloat32((i*vals.length+j)*4,vals[j],false);
  }
  return out;
}
