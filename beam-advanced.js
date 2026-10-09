// Simply supported beam: positive downward forces, positive clockwise applied couples.
// Internal positive bending is sagging. Deflection positive downward. x,L in m.
export function solveAdvancedBeam({L,pointLoads=[],patchLoads=[],moments=[],E=205000,I=80000000},steps=600){
 const valid=n=>Number.isFinite(n);
 if(![L,E,I].every(valid)||L<=0||E<=0||I<=0||!Number.isInteger(steps)||steps<10||steps>10000)throw Error('Invalid beam parameters');
 if(!pointLoads.every(q=>valid(q.P)&&q.P>=0&&valid(q.a)&&q.a>=0&&q.a<=L))throw Error('Invalid point load');
 if(!patchLoads.every(q=>valid(q.w)&&q.w>=0&&valid(q.a)&&valid(q.b)&&q.a>=0&&q.a<q.b&&q.b<=L))throw Error('Invalid partial UDL');
 if(!moments.every(q=>valid(q.M)&&valid(q.a)&&q.a>=0&&q.a<=L))throw Error('Invalid applied moment');
 const positive=x=>Math.max(0,x),H=x=>x>=0?1:0;
 const totalP=pointLoads.reduce((s,q)=>s+q.P,0)+patchLoads.reduce((s,q)=>s+q.w*(q.b-q.a),0);
 const RB=(pointLoads.reduce((s,q)=>s+q.P*q.a,0)+patchLoads.reduce((s,q)=>s+q.w*(q.b-q.a)*(q.a+q.b)/2,0)+moments.reduce((s,q)=>s+q.M,0))/L;
 const RA=totalP-RB,EI=E*I/1e9;
 function bending(x,side='right'){return RA*x-pointLoads.reduce((s,q)=>s+q.P*positive(x-q.a),0)-patchLoads.reduce((s,q)=>s+q.w*(positive(x-q.a)**2-positive(x-q.b)**2)/2,0)+moments.reduce((s,q)=>s+q.M*(side==='right'?x>=q.a:x>q.a),0)}
 function shear(x,side='right'){return RA-pointLoads.reduce((s,q)=>s+q.P*(side==='right'?x>=q.a:x>q.a),0)-patchLoads.reduce((s,q)=>s+q.w*(positive(x-q.a)-positive(x-q.b)),0)}
 function primitive(x){return RA*x**3/6-pointLoads.reduce((s,q)=>s+q.P*positive(x-q.a)**3/6,0)-patchLoads.reduce((s,q)=>s+q.w*(positive(x-q.a)**4-positive(x-q.b)**4)/24,0)+moments.reduce((s,q)=>s+q.M*positive(x-q.a)**2/2,0)}
 function slopePrimitive(x){return RA*x*x/2-pointLoads.reduce((s,q)=>s+q.P*positive(x-q.a)**2/2,0)-patchLoads.reduce((s,q)=>s+q.w*(positive(x-q.a)**3-positive(x-q.b)**3)/6,0)+moments.reduce((s,q)=>s+q.M*positive(x-q.a),0)}
 const C=primitive(L)/L;
 const at=(x,side='right')=>({x,V:shear(x,side),M:bending(x,side),y:1000*(C*x-primitive(x))/EI,theta:(C-slopePrimitive(x))/EI});
 const xs=new Set([0,L,...pointLoads.map(q=>q.a),...patchLoads.flatMap(q=>[q.a,q.b]),...moments.map(q=>q.a)]);
 for(let i=0;i<=steps;i++)xs.add(L*i/steps);
 const points=[...xs].sort((a,b)=>a-b).map(x=>at(x));
 const jumps=[...new Set([...pointLoads.map(q=>q.a),...moments.map(q=>q.a)])].sort((a,b)=>a-b).map(x=>({x,left:at(x,'left'),right:at(x,'right'),momentJump:moments.filter(q=>q.a===x).reduce((s,q)=>s+q.M,0)}));
 const candidates=[...points,...jumps.flatMap(j=>[j.left,j.right])];
 const maxMoment=Math.max(...candidates.map(p=>Math.abs(p.M)));
 const maxShear=Math.max(...candidates.map(p=>Math.abs(p.V)));
 const worstDeflection=points.reduce((a,b)=>Math.abs(b.y)>Math.abs(a.y)?b:a);
 return {RA,RB,points,jumps,maxMoment,maxShear,maxDeflection:Math.abs(worstDeflection.y),deflectionX:worstDeflection.x,totalLoad:totalP};
}
