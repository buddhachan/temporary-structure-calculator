// Euler-Bernoulli beam finite-element solver; kN, m, E N/mm2, I mm4.
// Vertical displacement is positive upward; applied loads P,w positive downward.
// Support: simple, cantilever-left, fixed-fixed, fixed-pinned.
export function solveSupportedBeam({L,E=205000,I=80000000,pointLoads=[],patchLoads=[],moments=[],support='simple',elements=80}){
 if(!Number.isFinite(L)||L<=0||!Number.isFinite(E)||E<=0||!Number.isFinite(I)||I<=0||!Number.isInteger(elements)||elements<4||elements>300)throw Error('Invalid beam parameters');
 if(!['simple','cantilever-left','fixed-fixed','fixed-pinned'].includes(support))throw Error('Invalid support');
 for(const q of pointLoads)if(!Number.isFinite(q.P)||q.P<0||!Number.isFinite(q.a)||q.a<0||q.a>L)throw Error('Invalid point load');
 for(const q of patchLoads)if(!Number.isFinite(q.w)||q.w<0||!Number.isFinite(q.a)||!Number.isFinite(q.b)||q.a<0||q.a>=q.b||q.b>L)throw Error('Invalid patch load');
 for(const q of moments)if(!Number.isFinite(q.M)||!Number.isFinite(q.a)||q.a<0||q.a>L)throw Error('Invalid moment');
 const xs=[0,L,...pointLoads.map(q=>q.a),...moments.map(q=>q.a),...patchLoads.flatMap(q=>[q.a,q.b])];for(let i=1;i<elements;i++)xs.push(L*i/elements);
 const nodes=[...new Set(xs)].sort((a,b)=>a-b),n=nodes.length,dof=2*n,EI=E*I/1e9;
 const K=Array.from({length:dof},()=>new Float64Array(dof)),F=new Float64Array(dof);
 const local=(h)=>{const c=EI/h**3;return [[12,6*h,-12,6*h],[6*h,4*h*h,-6*h,2*h*h],[-12,-6*h,12,-6*h],[6*h,2*h*h,-6*h,4*h*h]].map(row=>row.map(v=>v*c))};
 for(let i=0;i<n-1;i++){const a=nodes[i],b=nodes[i+1],h=b-a,k=local(h),ix=[2*i,2*i+1,2*i+2,2*i+3];for(let u=0;u<4;u++)for(let v=0;v<4;v++)K[ix[u]][ix[v]]+=k[u][v];for(const q of patchLoads)if(a>=q.a-1e-10&&b<=q.b+1e-10){const f=[-q.w*h/2,-q.w*h*h/12,-q.w*h/2,q.w*h*h/12];for(let u=0;u<4;u++)F[ix[u]]+=f[u]}}
 for(const q of pointLoads){const i=nodes.findIndex(x=>Math.abs(x-q.a)<1e-9);F[2*i]-=q.P}
 for(const q of moments){const i=nodes.findIndex(x=>Math.abs(x-q.a)<1e-9);F[2*i+1]-=q.M}
 const fixed=support==='simple'?[0,dof-2]:support==='cantilever-left'?[0,1]:support==='fixed-fixed'?[0,1,dof-2,dof-1]:[0,1,dof-2];
 const free=Array.from({length:dof},(_,i)=>i).filter(i=>!fixed.includes(i));const A=free.map(i=>Float64Array.from([...free.map(j=>K[i][j]),F[i]]));const m=free.length;
 for(let j=0;j<m;j++){let p=j;for(let i=j+1;i<m;i++)if(Math.abs(A[i][j])>Math.abs(A[p][j]))p=i;if(Math.abs(A[p][j])<1e-15)throw Error('Singular stiffness matrix');[A[p],A[j]]=[A[j],A[p]];const pivot=A[j][j];for(let k=j;k<=m;k++)A[j][k]/=pivot;for(let i=j+1;i<m;i++){const t=A[i][j];if(t===0)continue;for(let k=j;k<=m;k++)A[i][k]-=t*A[j][k]}}
 const u=new Float64Array(dof);for(let i=m-1;i>=0;i--){let v=A[i][m];for(let j=i+1;j<m;j++)v-=A[i][j]*u[free[j]];u[free[i]]=v}
 const reactions=fixed.map(i=>({dof:i,value:K[i].reduce((s,k,j)=>s+k*u[j],0)-F[i]}));
 const points=nodes.map((x,i)=>({x,y:-1000*u[2*i],theta:-u[2*i+1]}));
 const maxDeflection=Math.max(...points.map(p=>Math.abs(p.y)));
 const RA=reactions.find(r=>r.dof===0)?.value??0,RB=reactions.find(r=>r.dof===dof-2)?.value??0;
 const leftMoment=reactions.find(r=>r.dof===1)?.value??0,rightMoment=reactions.find(r=>r.dof===dof-1)?.value??0;
 const pos=x=>Math.max(0,x);
 const evaluate=(x,side='right')=>{
  const active=a=>side==='right'?x>=a:x>a;
  const V=RA-pointLoads.reduce((v,q)=>v+(active(q.a)?-q.P:0),0)-patchLoads.reduce((v,q)=>v+q.w*(pos(x-q.a)-pos(x-q.b)),0);
  const M=-leftMoment+RA*x-pointLoads.reduce((v,q)=>v+q.P*pos(x-q.a),0)-patchLoads.reduce((v,q)=>v+q.w*(pos(x-q.a)**2-pos(x-q.b)**2)/2,0)+moments.reduce((v,q)=>v+(active(q.a)?q.M:0),0);
  return {x,V,M};
 };
 const diagrams=nodes.map(x=>evaluate(x));
 const jumps=[...new Set([...pointLoads.map(q=>q.a),...moments.map(q=>q.a)])].sort((a,b)=>a-b).map(x=>({x,left:evaluate(x,'left'),right:evaluate(x,'right')}));
 const extrema=[...diagrams,...jumps.flatMap(j=>[j.left,j.right])];
 const maxShear=Math.max(...extrema.map(p=>Math.abs(p.V))),maxMoment=Math.max(...extrema.map(p=>Math.abs(p.M)));
 return {support,points,diagrams,jumps,reactions,RA,RB,leftMoment,rightMoment,maxDeflection,maxShear,maxMoment};
}
