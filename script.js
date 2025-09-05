/* =========================================================
   FactoriX — UPE + Small-Factor Sweep
   - Anti-blocage: boucles coopératives (await pause())
   - Boutons: Run / Stop / Reset
   - Journal défilant
   ========================================================= */

const TSP = "\u202f";               // espace fine
const ln10 = Math.log(10);
const state = { running:false, abort:false };

// ---------- Utils ----------
function groupBI(n){
  const s = n.toString();
  let out="", cnt=0;
  for(let i=s.length-1;i>=0;--i){
    out = s[i] + out; cnt++;
    if(cnt%3===0 && i!==0) out = TSP + out;
  }
  return out;
}
function isDigitsOnly(s){ return /^[0-9]+$/.test(s); }
function pause(ms=0){ return new Promise(res=>setTimeout(res,ms)); }
function logEl(){ return document.getElementById('log'); }
function addLine(html, cls=""){
  const div = document.createElement('div');
  div.className = `line ${cls}`;
  div.innerHTML = html;
  const L = logEl();
  L.appendChild(div);
  L.scrollTop = L.scrollHeight;
}
function clearLog(){ logEl().innerHTML = ""; }

// ---------- Parse input ----------
function parseBigIntExpr(raw){
  const s = raw.trim().replace(/\s+/g,"");
  if(s.length===0) return { ok:false, err:"Entrée vide." };
  if(isDigitsOnly(s)){
    const bi = BigInt(s);
    return { ok:true, mode:"big", big:bi, approxLog10:s.length-1 };
  }
  const m = /^10\^([0-9]+)(?:([+\-])([0-9]+))?$/.exec(s);
  if(m){
    const k = BigInt(m[1]);
    const op = m[2];
    const c = m[3] ? BigInt(m[3]) : 0n;
    if(k > 2000n){ // limite de sécurité
      return { ok:true, mode:"logonly", big:null, approxLog10:Number(k), adjust:{op,c} };
    }
    const tenPow = powBI(10n,k);
    let bi = tenPow;
    if(op==="+") bi = tenPow + c;
    if(op==="-") bi = tenPow - c;
    return { ok:true, mode:"big", big:bi, approxLog10:Number(k) };
  }
  return { ok:false, err:"Format non supporté. Utiliser chiffres ou 10^k (+/- c)." };
}
function powBI(a,e){
  let r=1n, b=a, ex=e;
  while(ex>0n){ if(ex&1n) r*=b; b*=b; ex>>=1n; }
  return r;
}
function lnBigIntApprox(bi){
  const s = bi.toString();
  const L = s.length;
  const head = s.slice(0, Math.min(18,L));
  const mant = Number(head)/Math.pow(10, head.length);
  return Math.log(mant) + L*ln10;
}

// ---------- Sieve small primes ----------
function primesUpTo(P){
  const n = Math.max(2, P|0);
  const sieve = new Uint8Array(n+1);
  const out = [];
  for(let i=2;i*i<=n;i++){
    if(!sieve[i]) for(let j=i*i;j<=n;j+=i) sieve[j]=1;
  }
  for(let i=2;i<=n;i++) if(!sieve[i]) out.push(i);
  return out;
}

// ---------- Miller–Rabin ----------
function modPow(base,exp,mod){
  let b = base%mod, e=exp, r=1n;
  while(e>0n){ if(e&1n) r=(r*b)%mod; b=(b*b)%mod; e>>=1n; }
  return r;
}
function isProbablePrime(n, rounds){
  if(n<2n) return false;
  const small=[2n,3n,5n,7n,11n,13n,17n,19n,23n,29n,31n,37n];
  for(const p of small){ if(n===p) return true; if(n%p===0n) return n===p; }
  let d=n-1n, s=0n; while((d&1n)===0n){ d>>=1n; s++; }
  function trial(a){
    if(a%n===0n) return true;
    let x = modPow(a,d,n);
    if(x===1n || x===n-1n) return true;
    for(let r=1n;r<s;r++){ x=(x*x)%n; if(x===n-1n) return true; }
    return false;
  }
  const fixed=[2n,3n,5n,7n,11n,13n,17n];
  let cnt=0;
  for(const a of fixed){ if(!trial(a)) return false; cnt++; if(cnt>=rounds) return true; }
  // pseudo-rand bases
  let seed=123456789;
  const span = n-3n;
  while(cnt<rounds){
    seed = (1103515245*seed + 12345) & 0x7fffffff;
    let a = 2n + (BigInt(seed) % (span>0n? span:1n));
    if(a>=n) a = (a % (n-2n)) + 2n;
    if(!trial(a)) return false;
    cnt++;
  }
  return true;
}

// ---------- Admissibility ----------
function admissibleAroundSqrt(x, t, S){
  if(state.abort) return false;
  const a = x - BigInt(t);
  const b = x + BigInt(t);
  if(a<2n || b<2n) return false;
  if((a%2n===0n && a!==2n) || (b%2n===0n && b!==2n)) return false;
  for(const s of S){
    const sb = BigInt(s);
    if(a===sb || b===sb) continue;
    if(a%sb===0n || b%sb===0n) return false;
  }
  return true;
}
function* symmetricOffsets(T){
  yield 0;
  for(let d=1; d<=T; d++){ yield +d; yield -d; }
}

// ---------- Small-factor sweep ----------
function quickSmallFactor(N, bound){
  const smalls = primesUpTo(bound);
  for(const p of smalls){
    if(state.abort) break;
    const pb = BigInt(p);
    if(N%pb===0n) return pb;
  }
  return null;
}

// ---------- Core factorization ----------
async function factorize(N, B, T, K, prefilter, mrRounds){
  const tStart = performance.now();

  // 0) trivial
  if(N%2n===0n) return { ok:true, p:2n, q:N/2n, phase:"even", ms:(performance.now()-tStart) };

  // 1) small-factor sweep
  if(B>0){
    addLine(`⛏️ Balayage petits facteurs ≤ ${B.toLocaleString('fr-FR')}…`, "muted");
    await pause(0);
    const sf = quickSmallFactor(N, Math.max(3, B|0));
    if(state.abort) return { ok:false, abort:true };
    if(sf){
      const ms = (performance.now()-tStart);
      return { ok:true, p:sf, q:N/sf, phase:"small", ms };
    }
    addLine(`— aucun petit facteur trouvé (≤ ${B.toLocaleString('fr-FR')}).`, "muted");
  }

  // 2) window around sqrt(N)
  const lnN = lnBigIntApprox(N);
  const P = Math.max(5, Math.floor(1.6 * lnN));   // petit crible par défaut
  const S = primesUpTo(P);
  const x = isqrt(N);

  addLine(`🔎 Fenêtre autour de √N ≈ ${groupBI(x)} | T=${T.toLocaleString('fr-FR')} | K=${K}`, "muted");
  await pause(0);

  let checked=0;
  for(let ring=0; ring<=K; ring++){
    if(state.abort) return { ok:false, abort:true };
    const R = ring * T;
    addLine(`— Anneau ${ring} (|u| ≤ ${R.toLocaleString('fr-FR')})…`, "muted");
    let steps=0;

    for(const u of symmetricOffsets(R)){
      if(state.abort) return { ok:false, abort:true };
      if(ring>0 && Math.abs(u) < (ring-1)*T) continue; // nouveau périmètre seulement
      const a = x + BigInt(u);
      const b = (u!==0)? (x - BigInt(u)) : null;

      // admissibilité optionnelle
      if(prefilter && !admissibleAroundSqrt(x, u, S)) continue;

      // test a
      checked++;
      if(N%a===0n){
        const ms = (performance.now()-tStart);
        return { ok:true, p:a, q:N/a, phase:"window", ring, u, checked, ms };
      }
      // test b
      if(b!==null){
        checked++;
        if(N%b===0n){
          const ms = (performance.now()-tStart);
          return { ok:true, p:b, q:N/b, phase:"window", ring, u:-u, checked, ms };
        }
      }

      steps++;
      if(steps % 5000 === 0){ // respiration UI
        addLine(`… ${steps.toLocaleString('fr-FR')} candidats testés sur anneau ${ring}`, "muted");
        await pause(0);
      }
    }
  }
  const ms = (performance.now()-tStart);
  return { ok:false, msg:"Aucun facteur dans la portée définie.", checked, ms };
}

// ---------- sqrt(BigInt) ----------
function isqrt(n){
  if(n<0n) throw new Error("isqrt negative");
  if(n<2n) return n;
  let x0 = 1n;
  let bit = 1n << (BigInt(n.toString(2).length - 1));
  while(bit>0n){
    const y = x0 + bit;
    if(y*y <= n) x0 = y;
    bit >>= 1n;
  }
  return x0;
}

// ---------- Wiring UI ----------
function getVal(id){ return document.getElementById(id).value; }
function getNum(id){ return parseInt(getVal(id),10) || 0; }
function getChk(id){ return document.getElementById(id).checked; }

async function onRun(){
  if(state.running){ addLine("Déjà en cours…", "warn"); return; }
  state.running = true; state.abort = false;

  const rawN = getVal('n');
  const B = getNum('B');
  const T = getNum('T');
  const K = getNum('K');
  const mr = getNum('mr') || 12;
  const pref = getChk('pref');
  const append = getChk('append');

  if(!append) clearLog();

  const parsed = parseBigIntExpr(rawN);
  if(!parsed.ok){ addLine(`Erreur entrée : ${parsed.err}`, "err"); state.running=false; return; }
  if(parsed.mode==="logonly"){
    addLine(`Mode aperçu (10^k très grand) — pas de boucle lourde.`, "warn");
    state.running=false; return;
  }
  const N = parsed.big;
  addLine(`<b>N</b> = ${groupBI(N)} <span class="muted">(≈ 10^${parsed.approxLog10})</span>`);

  try{
    const res = await factorize(N, B, T, K, pref, mr);
    if(res.abort){ addLine("⛔ Calcul interrompu (Stop).", "warn"); }
    else if(res.ok){
      addLine(`✅ <b>Facteurs</b> : p = ${groupBI(res.p)} · q = ${groupBI(res.q)} <span class="muted">[${res.phase}]</span>`, "ok");
      addLine(`⏱️ ${res.ms.toFixed(1)} ms`, "muted");
    }else{
      addLine(`❌ ${res.msg || "Échec."}`, "err");
      if(typeof res.checked === "number"){
        addLine(`Candidats testés ≈ ${res.checked?.toLocaleString('fr-FR')}`, "muted");
      }
      addLine(`⏱️ ${res.ms.toFixed(1)} ms`, "muted");
    }
  }catch(e){
    addLine(`Erreur inattendue : ${e && e.message ? e.message : String(e)}`, "err");
  }finally{
    state.running=false;
  }
}
function onStop(){ if(!state.running){ addLine("Rien à arrêter.", "warn"); return; } state.abort = true; }
function onReset(){
  state.abort = true;
  state.running = false;
  (document.getElementById('n').value = "");
  document.getElementById('B').value = 200000;
  document.getElementById('T').value = 20000;
  document.getElementById('K').value = 40;
  document.getElementById('mr').value = 12;
  document.getElementById('pref').checked = true;
  document.getElementById('append').checked = true;
  clearLog();
  addLine("Réinitialisé.", "muted");
}

document.addEventListener('DOMContentLoaded', ()=>{
  document.getElementById('run').addEventListener('click', onRun);
  document.getElementById('stop').addEventListener('click', onStop);
  document.getElementById('reset').addEventListener('click', onReset);
});
