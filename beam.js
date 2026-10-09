// Pure beam solver. Units: m, kN, N/mm², mm⁴.
export function solveBeam({L,loads=[],w=0,E=205000,I=80000000},steps=600){
 if(![L,w,E,I].every(Number.isFinite)||L<=0||w<0||E<=0||I<=0||!Number.isInteger(steps)||steps<2)throw Error('Invalid beam parameters');
 if(!loads.every(q=>Number.isFinite(q.P)&&q.P>=0&&Number.isFinite(q.a)&&q.a>=0&&q.a<=L))throw Error('Invalid point load');
 const RB=w*L/2+loads.reduce((s,q)=>s+q.P*q.a/L,0);
 const RA=w*L/2+loads.reduce((s,q)=>s+q.P*(L-q.a)/L,0);
 const EI=E*I,T=L*1000;
 function at(x){
  const X=x*1000;
  const V=RA-w*x-loads.reduce((s,q)=>s+(x>=q.a?q.P:0),0);
  const M=RA*x-w*x*x/2-loads.reduce((s,q)=>s+q.P*Math.max(0,x-q.a),0);
  let y=w*X*(T**3-2*T*X*X+X**3)/(24*EI);
  for(const {P,a} of loads){const A=a*1000,B=T-A; y+=X<=A?(P*1000)*B*X*(T*T-B*B-X*X)/(6*T*EI):(P*1000)*A*(T-X)*(T*T-A*A-(T-X)**2)/(6*T*EI)}
  return {x,V,M,y};
 }
 const xs=new Set([0,L,...loads.map(q=>q.a)]);for(let i=0;i<=steps;i++)xs.add(L*i/steps);
 const points=[...xs].sort((a,b)=>a-b).map(at);
 const moments=points.reduce((best,p)=>p.M>best.M?p:best,points[0]);
 const deflection=points.reduce((best,p)=>Math.abs(p.y)>Math.abs(best.y)?p:best,points[0]);
 const shearPositions=[0,L,...loads.map(q=>q.a)];
 const maxShear=Math.max(...shearPositions.flatMap(x=>{const before=RA-w*x-loads.reduce((s,q)=>s+(q.a<x?q.P:0),0);const after=RA-w*x-loads.reduce((s,q)=>s+(q.a<=x?q.P:0),0);return [Math.abs(before),Math.abs(after)]}));
 return {RA,RB,maxShear,maxMoment:moments.M,maxDeflection:Math.abs(deflection.y),deflectionX:deflection.x,points};
}
