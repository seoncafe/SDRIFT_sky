// Surface-brightness CCD limit; keep in sync with sensitivity.py.
// Exposure is accumulated usable time at fixed geometry, not elapsed orbit time.
// Binning is done in software after readout: read variance = N*n_reads*RN^2.
export const sensitivityDefaults={f_sys:.001,exposure_hours:10,bin_side:5,snr:3,
 diameter_cm:30,pixel_scale_arcsec:1.96,read_noise_e:1,dark_current_e_s:.05,n_reads:1,throughput:null};
const bands={UV1:[200,400],UV2:[250,400],UV3:[300,400],VIS:[400,700]};
export function sensitivityParameters(band='VIS',options={}){
 if(!bands[band])throw Error('Unknown sensitivity band');
 for(const k of Object.keys(options))if(!(k in sensitivityDefaults))throw Error('Unknown sensitivity setting: '+k);
 const p={...sensitivityDefaults,...options};if(p.throughput===null)p.throughput=band==='VIS'?.9:.5;
 for(const [k,v] of Object.entries(p)){
  if(!Number.isFinite(v))throw Error(k+' must be finite');
  if(['f_sys','read_noise_e','dark_current_e_s'].includes(k)?v<0:v<=0)throw Error(k+' is out of range');
 }
 if(p.f_sys>1||p.throughput>1)throw Error('f_sys and throughput must be fractions within 0–1');
 for(const k of ['bin_side','n_reads'])if(!Number.isInteger(p[k]))throw Error(k+' must be an integer');
 return p;
}
export function sensitivityCalculator(band='VIS',options={}){
 const p=sensitivityParameters(band,options),[lo,hi]=bands[band];
 const rate=1e-17/6.62607015e-27*Math.log(hi/lo)*Math.PI*(p.diameter_cm/2)**2*(p.pixel_scale_arcsec/206264.80624709636)**2*p.throughput;
 const t=p.exposure_hours*3600,n=p.bin_side**2,dark=n*p.dark_current_e_s*t,read=n*p.n_reads*p.read_noise_e**2;
 return bg=>{
  if(bg<0||(!Number.isFinite(bg)&&!Number.isNaN(bg)))throw Error('Invalid background');
  const sky=bg*rate*n*t,sys=p.f_sys*sky;
  const source=.5*p.snr**2*(1+Math.sqrt(1+4*(sky+dark+read+sys**2)/p.snr**2));
  const flux=source/(rate*n*t),mag=-2.5*Math.log10(flux*1e6/206264.80624709636**2/3631);
  return {flux,mag,sky,source};
 };
}
