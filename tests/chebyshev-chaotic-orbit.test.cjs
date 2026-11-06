const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const html = fs.readFileSync(path.join(__dirname, '../dynamical-systems/projects/chebyshev-cubic-kneading-diagram.html'), 'utf8');
const source = html.slice(html.indexOf('function chaoticOrbitSearch('), html.indexOf('function customOrbit()'));
const context = vm.createContext({});
vm.runInContext(source, context);
const search = context.chaoticOrbitSearch;
const map = (x, u, v) => (u - v) / 2 * (4 * x * x * x - 3 * x) + (u + v) / 2;
function forwardExponent(seed, u, v) {
 let x = seed, sum = 0;
 for (let n = 0; n < 8192; n++) {
  if (n >= 4096) sum += Math.log(Math.abs((u - v) * (6 * x * x - 1.5)));
  x = map(x, u, v);
 }
 return sum / 4096;
}
assert.ok(forwardExponent(-.5, .9, -.9) < 0);
assert.ok(forwardExponent(.5, .9, -.9) < 0);
for (const [u, v] of [[.79173279, -.65347389], [.65347389, -.79173279], [-.65347389, .79173279], [.8, -.8], [.9, -.9], [1, -1], [1.2, -1.2], [-1, 1], [1.8, -.2]]) {
 const result = search({u, v});
 assert.ok(!result.error, result.error);
 assert.equal(result.values.length, 104097);
 assert.equal(result.returnLengths.length, 2);
 assert.equal(result.method, 'kneading-symbolic-inverse');
 assert.ok(result.minBlockLambda > 0);
 assert.ok(result.kneadingDepth >= 16 && result.symbolicStates <= 8192);
 assert.ok(result.maxWidth < 1e-10);
 let error = 0;
 for (let i = 0; i < result.values.length - 1; i++) {
  const x = result.values[i];
  assert.ok(Number.isFinite(x) && Math.abs(x) < result.bound);
  error = Math.max(error, Math.abs(map(x, u, v) - result.values[i + 1]));
 }
 assert.ok(error < 1e-10, `Map residual ${error} at (${u}, ${v})`);
 for (const start of [0, 100000]) {
  let sum = 0;
  for (let i = start; i < start + 4096; i++) sum += Math.log(Math.abs((u - v) * (6 * result.values[i] ** 2 - 1.5)));
  assert.ok(sum / 4096 > 0);
 }
 for (let period = 1; period <= 64; period++) {
  assert.ok(result.values.subarray(256, 4352).some((x, i) => Math.abs(x - result.values[256 + i + period]) > 1e-8));
 }
}
assert.match(search({u: .3, v: .3}).error, /constant map/);
for (const [u,v] of [[.65,-.65],[.5,-.5],[.2,-.2]]) assert.match(search({u,v}).error, /does not establish/);
assert.match(search({u:NaN,v:0}).error, /finite parameters/);
console.log('Passed: expanding orbits, attracting and escaping critical-orbit cases, both orientations, map residuals, discard limit, and failure messages.');

// The reported failure point really is in a period-six attracting window.
const u=.79173279,v=-.65347389;
for(const seed of [-.5,.5]){
 assert.ok(forwardExponent(seed,u,v)<-.8);
 let x=seed;for(let i=0;i<10000;i++)x=map(x,u,v);
 let y=x;for(let i=0;i<6;i++)y=map(y,u,v);
 assert.ok(Math.abs(x-y)<1e-13);
}
const r=search({u,v}),other=search({u,v,offset:19});
assert.notEqual(r.symbols,other.symbols);
assert.ok(r.kneadingDepth===16, 'Regression point succeeds without deeper spatial preimages');

// Independent signed-order check on the returned numerical trace. Each
// position must strictly resolve its critical-value bounds within the prefix.
for(const [a,b] of [[u,v],[-v,-u],[-.65347389,.79173279],[1,-1],[-1,1]]){
 const result=search({u:a,v:b}),sgn=Math.sign(a-b),N=result.kneadingDepth;
 const code=x=>x<-.5?0:x===-.5?1:x<.5?2:x===.5?3:4;
 function itinerary(seed){let x=seed,out=[];for(let i=0;i<N;i++){out.push(code(x));x=map(x,a,b);}return out;}
 const upper=itinerary(Math.max(a,b)),lower=itinerary(Math.min(a,b));
 function compare(start,kneading){let parity=1;for(let j=0;j<N;j++){
  const c=code(result.values[start+j]);if(c!==kneading[j])return parity*Math.sign(c-kneading[j]);
  parity*=c===2?-sgn:sgn;
 }return 0;}
 for(let i=0;i<result.values.length-N;i++){
  const c=code(result.values[i]);
  if(c===2||(sgn>0?c===0:c===4))assert.equal(compare(i+1,upper),-1);
  if(c===2||(sgn>0?c===4:c===0))assert.equal(compare(i+1,lower),1);
 }
 if(Math.abs(a)===1&&a===-b)assert.ok(Math.abs(result.minBlockLambda-Math.log(3))<.005);
}
console.log('Passed: shrimp regression, periodic critical attractor, repeat-search variation, strict symbolic admissibility, and Chebyshev analytic Lyapunov benchmark.');
