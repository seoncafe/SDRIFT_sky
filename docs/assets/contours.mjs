// Join marching-square segments by shared grid edges before placing labels.
export function contourPaths(values,w,h,level,area){
 const nodes=[],edges=[],keys=new Map();
 function node(key,x,y){
  if(!keys.has(key)){keys.set(key,nodes.length);nodes.push({point:[area.x+x*area.w/w,area.y+y*area.h/h],edges:[]});}
  return keys.get(key);
 }
 for(let y=0;y<h;y++)for(let x=0;x<w;x++){
  const v=[values[y*(w+1)+x],values[y*(w+1)+x+1],values[(y+1)*(w+1)+x+1],values[(y+1)*(w+1)+x]];
  if(!v.every(Number.isFinite))continue;
  const corners=[[x,y],[x+1,y],[x+1,y+1],[x,y+1]],ids=[];
  const edgeKeys=[`h:${x}:${y}`,`v:${x+1}:${y}`,`h:${x}:${y+1}`,`v:${x}:${y}`];
  for(let k=0;k<4;k++){
   const j=(k+1)%4;if((v[k]<level)===(v[j]<level))continue;
   const f=(level-v[k])/(v[j]-v[k]);
   ids.push(node(edgeKeys[k],corners[k][0]+f*(corners[j][0]-corners[k][0]),corners[k][1]+f*(corners[j][1]-corners[k][1])));
  }
  // Preserve the existing bilinear-center choice for ambiguous saddle cells.
  if(ids.length===4&&((v.reduce((s,n)=>s+n,0)/4<level)!==(v[0]<level)))ids.push(ids.shift());
  for(let k=0;k+1<ids.length;k+=2){const a=ids[k],b=ids[k+1],e=edges.length;edges.push([a,b]);nodes[a].edges.push(e);nodes[b].edges.push(e);}
 }
 const used=new Uint8Array(edges.length),paths=[];
 function walk(start,edge){
  const points=[nodes[start].point];let n=start,length=0;
  while(edge!==undefined&&!used[edge]){
   used[edge]=1;const [a,b]=edges[edge],next=a===n?b:a,p=nodes[next].point,last=points[points.length-1];
   length+=Math.hypot(p[0]-last[0],p[1]-last[1]);points.push(p);n=next;
   if(n===start)break;
   edge=nodes[n].edges.find(e=>!used[e]);
  }
  if(length>1e-8)paths.push({points,length,closed:n===start});
 }
 // Open curves end at masks or the projection boundary; never bridge a gap.
 for(let n=0;n<nodes.length;n++)if(nodes[n].edges.length!==2)for(const e of nodes[n].edges)if(!used[e])walk(n,e);
 for(let e=0;e<edges.length;e++)if(!used[e])walk(edges[e][0],e);
 return paths;
}
function along(path,fraction){
 const target=path.length*fraction;let traversed=0;
 for(let i=1;i<path.points.length;i++){
  const a=path.points[i-1],b=path.points[i],d=Math.hypot(b[0]-a[0],b[1]-a[1]);
  if(d>0&&traversed+d>=target){const f=(target-traversed)/d;return [a[0]+f*(b[0]-a[0]),a[1]+f*(b[1]-a[1])];}
  traversed+=d;
 }
 return path.points[path.points.length-1];
}
// One mandatory label per curve. Callouts keep small disconnected curves readable.
export function contourLabels(curves,area,labelWidth=28,labelHeight=14){
 const labels=[],occupied=new Map(),cellSize=32,pad=2;
 const bounds={x:area.x+labelWidth/2,y:area.y+labelHeight/2,right:area.x+area.w-labelWidth/2,bottom:area.y+area.h-labelHeight/2};
 function cells(rect){const out=[];for(let y=Math.floor((rect.y-pad)/cellSize);y<=Math.floor((rect.y+rect.h+pad)/cellSize);y++)for(let x=Math.floor((rect.x-pad)/cellSize);x<=Math.floor((rect.x+rect.w+pad)/cellSize);x++)out.push(`${x}:${y}`);return out;}
 function overlaps(rect){
  let cost=0;const seen=new Set();for(const key of cells(rect))for(const i of occupied.get(key)||[]){
   if(seen.has(i))continue;seen.add(i);const r=labels[i].rect;
   cost+=Math.max(0,Math.min(rect.x+rect.w+pad,r.x+r.w)-Math.max(rect.x-pad,r.x))*Math.max(0,Math.min(rect.y+rect.h+pad,r.y+r.h)-Math.max(rect.y-pad,r.y));
  }return cost;
 }
 function place(curve,anchors){
  let best=null,bestScore=Infinity;
  for(const radius of [0,18,30,46,66,94]){
   for(const anchor of anchors)for(let angle=0;angle<(radius?8:1);angle++){
    const t=angle*Math.PI/4,x=Math.max(bounds.x,Math.min(bounds.right,anchor[0]+radius*Math.cos(t))),y=Math.max(bounds.y,Math.min(bounds.bottom,anchor[1]+radius*Math.sin(t)));
    const rect={x:x-labelWidth/2,y:y-labelHeight/2,w:labelWidth,h:labelHeight},collision=overlaps(rect),distance=Math.hypot(x-anchor[0],y-anchor[1]),score=collision*1000+distance;
    if(score<bestScore){bestScore=score;best={level:curve.level,x,y,anchor,rect,collision,callout:distance>8};}
    if(collision===0&&distance<=radius+2)break;
   }
   if(best.collision===0)break;
  }
  const index=labels.length;labels.push(best);for(const key of cells(best.rect)){if(!occupied.has(key))occupied.set(key,[]);occupied.get(key).push(index);}
 }
 const sorted=curves.slice().sort((a,b)=>b.length-a.length);
 // Complete mandatory labeling before spending space on repeated long-curve labels.
 for(const curve of sorted)place(curve,[.5,.3,.7,.15,.85].map(f=>along(curve,f)));
 for(const curve of sorted)if(curve.length>480){const count=Math.min(8,Math.floor(curve.length/240));for(let i=0;i<count;i++){const f=(i+.5)/count;if(Math.abs(f-.5)>.16)place(curve,[along(curve,f)]);}}
 return labels;
}
