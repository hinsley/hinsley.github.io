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
for (const [u, v] of [[.8, -.8], [.9, -.9], [1.2, -1.2], [-1, 1], [1.8, -.2]]) {
 const result = search({u, v});
 assert.ok(!result.error, result.error);
 assert.equal(result.values.length, 104097);
 assert.equal(result.returnLengths.length, 2);
 assert.ok(result.rate > 0);
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
assert.match(search({u: .65, v: -.65}).error, /does not mean/);
console.log('Passed: expanding orbits, attracting and escaping critical-orbit cases, both orientations, map residuals, discard limit, and failure messages.');
