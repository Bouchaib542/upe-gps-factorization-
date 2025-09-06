// script.js — FactoriXUPE solide (≤ 2^64-1). BigInt partout.

// ---------- Utils BigInt ----------
function isqrt(n){
  if (n < 0n) throw new Error("sqrt négative");
  if (n < 2n) return n;
  let x = n, y = (x + 1n) >> 1n;
  while (y < x){ x = y; y = (x + n / x) >> 1n; }
  return x;
}
function gcd(a,b){ while (b){ const t=a%b; a=b; b=t; } return a; }
function modPow(a,e,m){
  let r = 1n, x = a % m, ee = e;
  while (ee > 0n){
    if (ee & 1n) r = (r * x) % m;
    x = (x * x) % m;
    ee >>= 1n;
  }
  return r;
}

// ---------- Miller–Rabin déterministe (< 2^64) ----------
function isPrime64(n){
  if (n < 2n) return false;
  const small = [2n,3n,5n,7n,11n,13n,17n,19n,23n,29n];
  for (const p of small){ if (n % p === 0n) return n === p; }
  let d = n - 1n, s = 0n;
  while ((d & 1n) === 0n){ d >>= 1n; s++; }
  const bases = [2n,3n,5n,7n,11n,13n,17n];
  for (const a of bases){
    if (a % n === 0n) continue;
    let x = modPow(a, d, n);
    if (x === 1n || x === n-1n) continue;
    let ok = false;
    for (let r=1n; r<s; r++){
      x = (x*x) % n;
      if (x === n-1n){ ok = true; break; }
    }
    if (!ok) return false;
  }
  return true;
}

// ---------- Fermat (facteurs proches de √n) ----------
function fermatFactor(n, maxSteps = 2_000_000){
  let a = isqrt(n);
  if (a*a < n) a += 1n;
  for (let i=0; i<maxSteps; i++){
    const b2 = a*a - n;
    const b  = isqrt(b2);
    if (b*b === b2){
      const p = a - b, q = a + b;
      if (p>1n && q>1n && p*q === n) return [p, q];
    }
    a += 1n;
  }
  return null;
}

// ---------- Pollard-Rho (Brent) ----------
function rhoBrent(n){
  if (n % 2n === 0n) return 2n;
  const rand = () => (BigInt(Math.floor(Math.random()*1e9)) % (n-3n)) + 2n;
  let y = rand(), c = rand(), m = 1n<<6n;
  let g = 1n, r = 1n, q = 1n, x = 0n, ys = 0n;
  while (g === 1n){
    x = y;
    for (let i=0n; i<r; i++) y = (y*y + c) % n;
    let k = 0n;
    while (k < r && g === 1n){
      ys = y;
      const lim = (m < (r-k) ? m : (r-k));
      for (let i=0n; i<lim; i++){
        y = (y*y + c) % n;
        const diff = x > y ? x - y : y - x;
        q = (q * (diff % n)) % n;
      }
      g = gcd(q, n);
      k += m;
    }
    r <<= 1n;
  }
  if (g === n){
    do{
      ys = (ys*ys + c) % n;
      const diff = x > ys ? x - ys : ys - x;
      g = gcd(diff, n);
    } while (g === 1n);
  }
  return g;
}

// ---------- Factorisation récursive ----------
function factor(n, out){
  if (n === 1n) return;
  if (isPrime64(n)){ out.push(n); return; }
  // petit crible
  for (const p of [2n,3n,5n,7n,11n,13n,17n,19n,23n,29n,31n,37n,41n,43n,47n]){
    if (n % p === 0n){ factor(p, out); factor(n/p, out); return; }
  }
  // Fermat rapide
  const ff = fermatFactor(n, 2_000_000);
  if (ff){ factor(ff[0], out); factor(ff[1], out); return; }
  // Pollard-Rho (Brent)
  let d = n;
  while (d === n) d = rhoBrent(n);
  factor(d, out); factor(n/d, out);
}

export function factorizeBigInt(N){
  const res = [];
  factor(N, res);
  return res;
}

// ---------- Hook UI (si présent dans index.html) ----------
function $(s){ return document.querySelector(s); }
function log(s){ const out = $('#out'); if (out) out.textContent += s + "\n"; }

if (typeof window !== 'undefined'){
  window.addEventListener('DOMContentLoaded', ()=>{
    const nin = $('#nin'), go = $('#go'), out = $('#out');
    if (!nin || !go || !out) return;
    go.onclick = ()=>{
      out.textContent = "";
      const txt = nin.value.trim();
      if (!/^\d+$/.test(txt)){ log("Entrée invalide. Saisir un entier décimal."); return; }
      const N = BigInt(txt);
      const t0 = performance.now();
      const fs = factorizeBigInt(N).sort((a,b)=> (a<b?-1:1));
      const t1 = performance.now();
      const prod = fs.reduce((p,x)=> p*x, 1n);
      log("N = " + N.toString());
      log("Facteurs ("+fs.length+"): " + fs.map(x=>x.toString()).join(" × "));
      log("Vérif: produit == N  →  " + (prod === N));
      log("Temps: " + (t1 - t0).toFixed(1) + " ms");
    };
  });
}
