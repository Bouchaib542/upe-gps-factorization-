/* =========================================================
   Unified Prime Equation (UPE) & GPS Factorization — script
   - Prime near X
   - Goldbach (E = 2x)
   - Factorization (N = p·q) with sieve-guided √N window
   ========================================================= */

/* -------------------- Utilities -------------------- */

// Thin space for pretty prints
const TSP = "\u202f";

function groupBI(n) {
  // BigInt -> grouped string
  const s = n.toString();
  let out = "";
  let cnt = 0;
  for (let i = s.length - 1; i >= 0; --i) {
    out = s[i] + out;
    cnt++;
    if (cnt % 3 === 0 && i !== 0) out = TSP + out;
  }
  return out;
}

function isDigitsOnly(s) {
  return /^[0-9]+$/.test(s);
}

// Parse "10^k", "10^k +/- c", or plain digits. Returns { ok, big, approxLog10, mode }
function parseBigIntExpr(raw) {
  const s = raw.trim().replace(/\s+/g, "");
  if (isDigitsOnly(s)) {
    const bi = BigInt(s);
    return { ok: true, big: bi, approxLog10: s.length - 1, mode: "big" };
  }
  // pattern: 10^k (+/- c) ?  (k, c are digits)
  const m = /^10\^([0-9]+)(?:([+\-])([0-9]+))?$/.exec(s);
  if (m) {
    const k = BigInt(m[1]);
    const op = m[2];
    const c = m[3] ? BigInt(m[3]) : 0n;
    if (k > 20000n) {
      // Too huge to build exact BigInt; preview/log mode only
      return { ok: true, big: null, approxLog10: Number(k), mode: "logonly", adjust: { op, c } };
    }
    // Build 10^k exactly
    let tenPow = 1n;
    const ten = 10n;
    let kk = k;
    // exponentiation by squaring
    function powBI(a, e) {
      let r = 1n;
      let base = a;
      let exp = e;
      while (exp > 0n) {
        if (exp & 1n) r *= base;
        base *= base;
        exp >>= 1n;
      }
      return r;
    }
    tenPow = powBI(ten, kk);
    let bi = tenPow;
    if (op === "+") bi = tenPow + c;
    if (op === "-") bi = tenPow - c;
    return { ok: true, big: bi, approxLog10: Number(k), mode: "big" };
  }
  return { ok: false, err: "Unsupported input. Use digits or 10^k (+/- c)." };
}

// Approx ln(BigInt) using leading digits (string trick)
function lnBigIntApprox(bi) {
  const s = bi.toString();
  const L = s.length;
  // take first 18 digits to form a float mantissa
  const head = s.slice(0, Math.min(18, L));
  const mant = Number(head) / Math.pow(10, head.length);
  const ln10 = Math.log(10);
  return Math.log(mant) + (L * ln10);
}

// floor sqrt for BigInt
function isqrt(n) {
  if (n < 0n) throw new Error("isqrt of negative");
  if (n < 2n) return n;
  // Newton
  let x0 = 1n << (BigInt(n.toString().length) * 2n); // rough
  let x1 = (x0 + n / x0) >> 1n;
  while (x1 < x0) {
    x0 = x1;
    x1 = (x0 + n / x0) >> 1n;
  }
  return x0;
}

/* -------------------- Small primes up to P -------------------- */

function primesUpTo(P) {
  const n = Math.max(2, P|0);
  const sieve = new Uint8Array(n + 1);
  const out = [];
  for (let i = 2; i * i <= n; i++) {
    if (!sieve[i]) {
      for (let j = i * i; j <= n; j += i) sieve[j] = 1;
    }
  }
  for (let i = 2; i <= n; i++) if (!sieve[i]) out.push(i);
  return out;
}

/* -------------------- Miller–Rabin (BigInt) -------------------- */

function modPow(base, exp, mod) {
  let b = base % mod;
  let e = exp;
  let r = 1n;
  while (e > 0n) {
    if (e & 1n) r = (r * b) % mod;
    b = (b * b) % mod;
    e >>= 1n;
  }
  return r;
}

function isProbablePrime(n, rounds) {
  if (n < 2n) return false;
  const small = [2n,3n,5n,7n,11n,13n,17n,19n,23n,29n,31n,37n];
  for (const p of small) {
    if (n === p) return true;
    if (n % p === 0n) return n === p;
  }
  // write n-1 = d*2^s
  let d = n - 1n, s = 0n;
  while ((d & 1n) === 0n) { d >>= 1n; s++; }
  function trial(a) {
    if (a % n === 0n) return true;
    let x = modPow(a, d, n);
    if (x === 1n || x === n - 1n) return true;
    for (let r = 1n; r < s; r++) {
      x = (x * x) % n;
      if (x === n - 1n) return true;
    }
    return false;
  }
  // use fixed bases first, then a few pseudo-random
  const fixed = [2n,3n,5n,7n,11n,13n,17n];
  let cnt = 0;
  for (const a of fixed) {
    if (!trial(a)) return false;
    cnt++; if (cnt >= rounds) return true;
  }
  // pseudo-random (deterministic seed)
  let seed = 123456789;
  function rnd() {
    seed = (1103515245 * seed + 12345) & 0x7fffffff;
    return BigInt(seed);
  }
  while (cnt < rounds) {
    const a = 2n + (rnd() % (n - 3n));
    if (!trial(a)) return false;
    cnt++;
  }
  return true;
}

/* -------------------- Admissibility (sieve) -------------------- */

function admissibleSingle(X, u, S) {
  const Xu = X + BigInt(u);
  if (Xu < 2n) return false;
  if (Xu % 2n === 0n) return Xu === 2n;
  for (const s of S) {
    const sb = BigInt(s);
    if (Xu === sb) return true;
    if (Xu % sb === 0n) return false;
  }
  return true;
}

function admissibleGoldbach(x, t, S) {
  const a = x - BigInt(t);
  const b = x + BigInt(t);
  if (a < 2n || b < 2n) return false;
  if ((a % 2n === 0n && a !== 2n) || (b % 2n === 0n && b !== 2n)) return false;
  for (const s of S) {
    const sb = BigInt(s);
    if (a === sb || b === sb) continue;
    if (a % sb === 0n || b % sb === 0n) return false;
  }
  return true;
}

function* symmetricOffsets(T) {
  yield 0;
  for (let d = 1; d <= T; d++) {
    yield +d;
    yield -d;
  }
}

/* -------------------- Prime near X -------------------- */

function runPrimeNearX(rawX, c1, c2, Poverride, Toverride, mrRounds, outEl) {
  const parsed = parseBigIntExpr(rawX);
  if (!parsed.ok) {
    outEl.textContent = "Error: " + parsed.err;
    return;
  }
  if (parsed.mode === "logonly") {
    const ln10 = Math.log(10);
    const lnN = parsed.approxLog10 * ln10;
    const P = Math.max(5, Math.floor(c1 * lnN));
    const T = Math.max(8, Math.floor(c2 * (lnN * lnN)));
    outEl.textContent =
      `Log Mode Preview\n` +
      `log₁₀ N ≈ ${parsed.approxLog10}\n` +
      `ln N ≈ ${lnN.toFixed(6)}\n` +
      `Small-prime cutoff P ≈ ${P}\n` +
      `Window radius |u| ≤ T ≈ ${T}\n` +
      `Recommendation: use external verifier for primality at this scale.`;
    return;
  }
  const X = parsed.big;
  const lnN = lnBigIntApprox(X);
  const P = Poverride ? Math.max(3, Poverride|0) : Math.max(5, Math.floor(c1 * lnN));
  const T = Toverride ? Math.max(2, Toverride|0) : Math.max(8, Math.floor(c2 * (lnN * lnN)));
  const S = primesUpTo(P);

  const t0 = performance.now();
  let checked = 0;
  for (const u of symmetricOffsets(T)) {
    if (!admissibleSingle(X, u, S)) continue;
    checked++;
    const cand = X + BigInt(u);
    const prime = isProbablePrime(cand, mrRounds);
    if (prime) {
      const ms = (performance.now() - t0).toFixed(2);
      outEl.textContent =
        `Prime near X\n` +
        `Input X = ${groupBI(X)}\n` +
        `P = ${P}, T = ${T}, MR rounds = ${mrRounds}\n` +
        `First admissible found at u = ${u}\n` +
        `Prime = ${groupBI(cand)}\n` +
        `Admissibles tested = ${checked} (Δ_step = ${checked-1})\n` +
        `Time = ${ms} ms`;
      return;
    }
    if (checked >= 1 + 2) break; // bounded correction (≤ 2) as per UPE heuristic
  }
  const ms = (performance.now() - t0).toFixed(2);
  outEl.textContent =
    `No prime found within bounded correction.\n` +
    `Consider increasing T or MR rounds.\n` +
    `Checked = ${checked}, Time = ${ms} ms`;
}

/* -------------------- Goldbach (E = 2x) -------------------- */

function runGoldbach(rawE, c1, c2, Poverride, Toverride, mrRounds, outEl) {
  const parsed = parseBigIntExpr(rawE);
  if (!parsed.ok) { outEl.textContent = "Error: " + parsed.err; return; }

  if (parsed.mode === "logonly") {
    const ln10 = Math.log(10);
    const lnN = parsed.approxLog10 * ln10;
    const P = Math.max(5, Math.floor(c1 * lnN));
    const T = Math.max(8, Math.floor(c2 * (lnN * lnN)));
    outEl.textContent =
      `Log Mode Preview\n` +
      `log₁₀ E ≈ ${parsed.approxLog10}\n` +
      `ln E ≈ ${lnN.toFixed(6)}\n` +
      `Small-prime cutoff P ≈ ${P}\n` +
      `Window radius |t| ≤ T ≈ ${T}\n` +
      `Recommendation: use external verifier for primality at this scale.`;
    return;
  }
  const E = parsed.big;
  if (E % 2n !== 0n) { outEl.textContent = "E must be even."; return; }

  const x = E / 2n;
  const lnE = lnBigIntApprox(E);
  const P = Poverride ? Math.max(3, Poverride|0) : Math.max(5, Math.floor(c1 * lnE));
  const T = Toverride ? Math.max(2, Toverride|0) : Math.max(8, Math.floor(c2 * (lnE * lnE)));
  const S = primesUpTo(P);

  const t0 = performance.now();
  let checked = 0;
  for (const t of symmetricOffsets(T)) {
    if (!admissibleGoldbach(x, t, S)) continue;
    checked++;
    const a = x - BigInt(t);
    const b = x + BigInt(t);
    if (isProbablePrime(a, mrRounds) && isProbablePrime(b, mrRounds)) {
      const ms = (performance.now() - t0).toFixed(2);
      outEl.textContent =
        `Goldbach Pair\n` +
        `Input E = ${groupBI(E)}, x = E/2 = ${groupBI(x)}\n` +
        `P = ${P}, T = ${T}, MR rounds = ${mrRounds}\n` +
        `First admissible t = ${t}\n` +
        `Pair: (${groupBI(a)}, ${groupBI(b)})\n` +
        `Admissibles tested = ${checked} (Δ_step = ${checked-1})\n` +
        `Time = ${ms} ms`;
      return;
    }
    if (checked >= 1 + 2) break; // bounded correction (≤ 2)
  }
  const ms = (performance.now() - t0).toFixed(2);
  outEl.textContent =
    `No pair found within bounded correction.\n` +
    `Consider increasing T or MR rounds.\n` +
    `Checked = ${checked}, Time = ${ms} ms`;
}

/* -------------------- Factorization (N = p·q) -------------------- */

function runFactorization(rawN, c1, c2, K, usePrefilter, mrRounds, outEl) {
  const parsed = parseBigIntExpr(rawN);
  if (!parsed.ok) { outEl.textContent = "Error: " + parsed.err; return; }

  if (parsed.mode === "logonly") {
    const ln10 = Math.log(10);
    const lnN = parsed.approxLog10 * ln10;
    const P = Math.max(5, Math.floor(c1 * lnN));
    const T = Math.max(8, Math.floor(c2 * (lnN * lnN)));
    outEl.textContent =
      `Log Mode Preview\n` +
      `log₁₀ N ≈ ${parsed.approxLog10}\n` +
      `ln N ≈ ${lnN.toFixed(6)}\n` +
      `Small-prime cutoff P ≈ ${P}\n` +
      `Window radius |u| ≤ T ≈ ${T}\n` +
      `Suggestion: if wide gaps suspected, use QS-on-tracks offline.`;
    return;
  }

  const N = parsed.big;
  if (N < 4n) { outEl.textContent = "N must be ≥ 4."; return; }
  const lnN = lnBigIntApprox(N);
  const P = Math.max(5, Math.floor(c1 * lnN));
  const T = Math.max(8, Math.floor(c2 * (lnN * lnN)));
  const S = primesUpTo(P);

  const X = isqrt(N); // anchor near √N

  // Optional prefilter: require both (X±u) coprime to all primes ≤ P (same as “admissibleSingle”)
  function passesPrefilter(y) {
    if (!usePrefilter) return true;
    if (y < 2n) return false;
    for (const s of S) {
      const sb = BigInt(s);
      if (y === sb) return true;
      if (y % sb === 0n) return false;
    }
    return true;
  }

  const t0 = performance.now();
  let checked = 0;
  // Explore ring by ring: R = ring*T, for ring=0..K
  for (let ring = 0; ring <= K; ring++) {
    const R = ring * T;
    const from = -R, to = +R;
    for (const u of symmetricOffsets(R)) {
      if (Math.abs(u) < (ring-1)*T) continue; // only the new ring
      const candA = X + BigInt(u);
      const candB = (u !== 0) ? (X - BigInt(u)) : null;
      if (passesPrefilter(candA)) {
        checked++;
        if (N % candA === 0n) {
          const p = candA;
          const q = N / candA;
          const ms = (performance.now() - t0).toFixed(2);
          outEl.textContent =
            `Factorization\n` +
            `N = ${groupBI(N)}\n` +
            `√N ≈ ${groupBI(X)} | Window T = ${T}, Rings K = ${K}, P = ${P}\n` +
            `Hit at u = ${u} → factor = ${groupBI(p)}, partner = ${groupBI(q)}\n` +
            `Admissible tests = ${checked}\n` +
            `Time = ${ms} ms`;
          return;
        }
      }
      if (candB !== null && passesPrefilter(candB)) {
        checked++;
        if (N % candB === 0n) {
          const p = candB;
          const q = N / candB;
          const ms = (performance.now() - t0).toFixed(2);
          outEl.textContent =
            `Factorization\n` +
            `N = ${groupBI(N)}\n` +
            `√N ≈ ${groupBI(X)} | Window T = ${T}, Rings K = ${K}, P = ${P}\n` +
            `Hit at u = ${-u} → factor = ${groupBI(p)}, partner = ${groupBI(q)}\n` +
            `Admissible tests = ${checked}\n` +
            `Time = ${ms} ms`;
          return;
        }
      }

      // Bounded correction shortcut for small rings
      if (ring === 0 && checked >= 3) {
        const ms = (performance.now() - t0).toFixed(2);
        outEl.textContent =
          `No factor within bounded correction at ring 0.\n` +
          `Consider increasing Rings (K) or using QS-on-tracks offline.\n` +
          `Checked = ${checked}, Time = ${ms} ms`;
        return;
      }
    }
  }
  const ms = (performance.now() - t0).toFixed(2);
  outEl.textContent =
    `No factor found within K rings.\n` +
    `Try larger K, or switch to QS-on-tracks for wide gaps.\n` +
    `Checked = ${checked}, Time = ${ms} ms`;
}

/* -------------------- Tabs & Wiring -------------------- */

function setupTabs() {
  const tabs = Array.from(document.querySelectorAll(".tab"));
  const panels = {
    prime: document.querySelector("#panel-prime"),
    goldbach: document.querySelector("#panel-goldbach"),
    factor: document.querySelector("#panel-factor"),
  };
  tabs.forEach(btn => {
    btn.addEventListener("click", () => {
      tabs.forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
      const key = btn.dataset.tab;
      Object.values(panels).forEach(p => p.classList.remove("active"));
      panels[key].classList.add("active");
      if (location.hash !== `#${key}`) history.replaceState(null, "", `#${key}`);
    });
  });
  // Deep link tab
  if (location.hash) {
    const key = location.hash.replace("#", "");
    const target = document.querySelector(`.tab[data-tab="${key}"]`);
    if (target) target.click();
  }
}

function setupPrime() {
  const xEl = document.querySelector("#prime-x");
  const c1El = document.querySelector("#prime-c1");
  const c2El = document.querySelector("#prime-c2");
  const PEl  = document.querySelector("#prime-P-override");
  const TEl  = document.querySelector("#prime-T-override");
  const mrEl = document.querySelector("#prime-mr");
  const run  = document.querySelector("#prime-run");
  const clr  = document.querySelector("#prime-clear");
  const out  = document.querySelector("#prime-out");

  run.addEventListener("click", () => {
    out.textContent = "Searching prime near X...";
    runPrimeNearX(
      xEl.value,
      parseFloat(c1El.value || "1.6"),
      parseFloat(c2El.value || "2.0"),
      PEl.value ? parseInt(PEl.value, 10) : 0,
      TEl.value ? parseInt(TEl.value, 10) : 0,
      parseInt(mrEl.value, 10) || 12,
      out
    );
  });
  clr.addEventListener("click", () => out.textContent = "");
}

function setupGoldbach() {
  const eEl = document.querySelector("#gb-e");
  const c1El = document.querySelector("#gb-c1");
  const c2El = document.querySelector("#gb-c2");
  const PEl  = document.querySelector("#gb-P-override");
  const TEl  = document.querySelector("#gb-T-override");
  const mrEl = document.querySelector("#gb-mr");
  const run  = document.querySelector("#gb-run");
  const clr  = document.querySelector("#gb-clear");
  const out  = document.querySelector("#gb-out");

  run.addEventListener("click", () => {
    out.textContent = "Scanning around x = E/2...";
    runGoldbach(
      eEl.value,
      parseFloat(c1El.value || "1.6"),
      parseFloat(c2El.value || "2.0"),
      PEl.value ? parseInt(PEl.value, 10) : 0,
      TEl.value ? parseInt(TEl.value, 10) : 0,
      parseInt(mrEl.value, 10) || 12,
      out
    );
  });
  clr.addEventListener("click", () => out.textContent = "");
}

function setupFactor() {
  const nEl = document.querySelector("#fac-n");
  const c1El = document.querySelector("#fac-c1");
  const c2El = document.querySelector("#fac-c2");
  const kEl  = document.querySelector("#fac-k");
  const pref = document.querySelector("#fac-pref");
  const mrEl = document.querySelector("#fac-mr");
  const run  = document.querySelector("#fac-run");
  const clr  = document.querySelector("#fac-clear");
  const out  = document.querySelector("#fac-out");

  run.addEventListener("click", () => {
    out.textContent = "Factoring near √N...";
    runFactorization(
      nEl.value,
      parseFloat(c1El.value || "1.6"),
      parseFloat(c2El.value || "2.0"),
      parseInt(kEl.value, 10) || 40,
      (pref.value || "on") === "on",
      parseInt(mrEl.value, 10) || 12,
      out
    );
  });
  clr.addEventListener("click", () => out.textContent = "");
}

document.addEventListener("DOMContentLoaded", () => {
  setupTabs();
  setupPrime();
  setupGoldbach();
  setupFactor();
});
```0
