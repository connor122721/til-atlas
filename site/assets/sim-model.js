/* Stochastic model of infused TIL clones after adoptive transfer.
 *
 * Each infused clone i carries three phenotype compartments:
 *   S  stem-like / memory-progenitor   (TCF7+, CD39-CD69-  — Krishna et al. 2020 "DN")
 *   E  effector                        (GZMB+, kills tumor)
 *   X  terminally exhausted            (CD39+CD69+, TOX+   — Krishna "DP"; kills weakly)
 * Clones are tumor-reactive (r_i = 1, functional avidity a_i) or bystanders (r_i = 0).
 *
 * Per day, with crowding c = max(0, 1 - N / Ncap):
 *   antigen drive      A_i(t) = r_i a_i * Tu / (Tu + K_T)
 *   homeostatic drive  H(t)   = h0 exp(-t / tau_h)            (space after lymphodepletion: IL-7/IL-15)
 *   IL-2 drive         I(t)   = il2 for t < il2Days
 *   S: self-renew  rho_S (H + I + sA A) c          (S -> S + S)
 *      asymmetric  rho_A A                         (S -> S + E : the stem pool is an effector factory)
 *      lose stemness delta A,  die mu_S
 *   E: divide  rho_E (H + I + A) c,   E->X  eps A,   die mu_E   (effectors are short-lived)
 *   X: die mu_X                                                 (exhausted cells barely divide)
 * Tumor (relative burden Tu, 1 = baseline):
 *   dTu/dt = g Tu (1 - Tu / Tmax) - kmax Tu * K / (K + Kh),   K = sum_i r_i a_i (E_i + 0.3 X_i)
 * Endogenous T cells regrow logistically after lymphodepletion and dilute the infused clones in blood.
 *
 * Cells are integers in units of 10^6 and every event count is a Poisson draw (tau-leaping, dt = 0.25 d),
 * so small clones can go extinct by chance — that's what makes clonal dynamics interesting.
 * Parameters are illustrative, chosen to reproduce qualitative patterns from the literature, not fitted.
 */
(function (root) {
  function rng(seed) { // mulberry32
    let a = seed >>> 0;
    return () => { a = (a + 0x6d2b79f5) >>> 0; let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
  }
  function normal(R) { let u = 0, v = 0; while (!u) u = R(); while (!v) v = R(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
  function poisson(R, lam) {
    if (lam <= 0) return 0;
    if (lam < 30) { const L = Math.exp(-lam); let k = 0, p = 1; do { k++; p *= R(); } while (p > L); return k - 1; }
    return Math.max(0, Math.round(lam + Math.sqrt(lam) * normal(R)));
  }

  const DEFAULTS = {
    seed: 7, days: 365, dt: 0.25,
    nClones: 80,          // clones tracked individually
    productCells: 3e4,    // x10^6 cells = 3e10
    clonality: 1.1,       // Zipf exponent for clone sizes in the product
    stemFrac: 0.20,       // fraction of each clone's cells that are stem-like (S)
    reactiveFrac: 0.25,   // fraction of clones that recognize the tumor
    avidity: 1.0,         // mean functional avidity of reactive clones
    lymphodep: true, il2: true, il2Days: 3,
    tumorGrowth: 0.035,   // per day (doubling ~ 20 d)
    resistance: 1.0,      // multiplies Kh: higher = harder to kill (antigen loss, MHC-I loss, TME)
  };
  const PARAMS = { rhoS: 0.5, sA: 0.05, rhoA: 0.8, delta: 0.03, rhoE: 0.3, eps: 0.25, muS: 0.015, muE: 0.15, muX: 0.07,
    KT: 0.05, h0: 0.5, tauH: 7, il2: 0.45, Ncap: 1.2e5, Tmax: 6, kmax: 0.35, Kh: 12000,
    endo0: 400, endoR: 0.07, endoK: 3e4, cure: 1e-4 };

  function simulate(opts = {}) {
    const o = { ...DEFAULTS, ...opts }, R = rng(o.seed), P = { ...PARAMS, ...(opts.params || {}) };
    // ---- build the infusion product
    let clones;
    if (o.repertoire) {
      // REAL INPUT: clone sizes and S/E/X composition from a patient's tumor CD8 repertoire, scaled to the
      // infused cell number. Avidity is not measured, so predicted-reactive clones draw it around o.avidity.
      const tot = o.repertoire.reduce((a, c) => a + c.S + c.E + c.X, 0);
      const k = o.productCells / tot;
      clones = o.repertoire.map((c, i) => {
        const S = Math.round(c.S * k), E = Math.round(c.E * k), X = Math.max(Math.round(c.X * k), S + E ? 0 : 1);
        const a = c.reactive ? o.avidity * Math.exp(0.45 * normal(R)) : 0;
        return { id: i + 1, label: c.id, reactive: c.reactive, avidity: a, cls: c.reactive ? "high" : "bystander",
          S, E, X, infused: S + E + X, stem0: (c.S) / (c.S + c.E + c.X), inBlood: c.inBlood, inNAT: c.inNAT };
      });
    } else {
      const w = Array.from({ length: o.nClones }, (_, k) => Math.pow(k + 1, -o.clonality));
      const wsum = w.reduce((a, b) => a + b, 0);
      const order = Array.from({ length: o.nClones }, (_, k) => k).sort(() => R() - 0.5); // which ranks are reactive
      const nReactive = Math.round(o.reactiveFrac * o.nClones);
      clones = w.map((wi, k) => {
        const n = Math.max(1, Math.round(o.productCells * wi / wsum));
        const reactive = order.indexOf(k) < nReactive;
        const a = reactive ? o.avidity * Math.exp(0.45 * normal(R)) : 0;
        const s = Math.min(0.95, Math.max(0, o.stemFrac * Math.exp(0.4 * normal(R))));
        const S = Math.round(n * s), X = Math.round((n - S) * 0.55), E = n - S - X;
        return { id: k + 1, reactive, avidity: a, cls: !reactive ? "bystander" : a >= o.avidity ? "high" : "low",
          S, E, X, infused: n, stem0: s };
      });
    }

    // ---- integrate
    const steps = Math.round(o.days / o.dt), every = Math.round(1 / o.dt);
    let Tu = 1, endo = o.lymphodep ? P.endo0 : P.endoK * 0.8;
    const tl = { t: [], tumor: [], S: [], E: [], X: [], endo: [], clones: clones.map(() => []) };
    const record = (t) => {
      let S = 0, E = 0, X = 0;
      clones.forEach((c, i) => { S += c.S; E += c.E; X += c.X; tl.clones[i].push(c.S + c.E + c.X); });
      tl.t.push(t); tl.tumor.push(Tu); tl.S.push(S); tl.E.push(E); tl.X.push(X); tl.endo.push(endo);
    };
    record(0);
    for (let step = 1; step <= steps; step++) {
      const t = step * o.dt, dt = o.dt;
      const H = (o.lymphodep ? P.h0 : 0.05) * Math.exp(-t / P.tauH);
      const I = o.il2 && t < o.il2Days ? P.il2 : 0;
      let N = endo, K = 0;
      clones.forEach((c) => { N += c.S + c.E + c.X; K += c.avidity * (c.E + 0.3 * c.X); });
      const crowd = Math.max(0, 1 - N / P.Ncap), tumorSat = Tu / (Tu + P.KT);
      for (const c of clones) {
        const A = c.avidity * tumorSat;
        const bS = poisson(R, P.rhoS * (H + I + P.sA * A) * crowd * c.S * dt);
        const gE = poisson(R, P.rhoA * A * c.S * dt);                   // asymmetric: new effectors, S unchanged
        const dS = Math.min(c.S, poisson(R, P.muS * c.S * dt));
        const tSE = Math.min(c.S - dS, poisson(R, P.delta * A * c.S * dt));
        const bE = poisson(R, P.rhoE * (H + I + A) * crowd * c.E * dt);
        const dE = Math.min(c.E, poisson(R, P.muE * c.E * dt));
        const tEX = Math.min(c.E - dE, poisson(R, P.eps * A * c.E * dt));
        const dX = Math.min(c.X, poisson(R, P.muX * c.X * dt));
        c.S += bS - dS - tSE; c.E += gE + bE - dE - tEX + tSE; c.X += tEX - dX;
      }
      endo += P.endoR * endo * (1 - endo / P.endoK) * dt;
      Tu += (o.tumorGrowth * Tu * (1 - Tu / P.Tmax) - P.kmax * Tu * K / (K + P.Kh * o.resistance)) * dt;
      if (Tu < P.cure) Tu = 0; // below ~1 lesion's worth of cells: cured, and antigen drive stops
      if (step % every === 0) record(t);
    }
    return { opts: o, clones, timeline: tl, readouts: readouts(tl, clones, R) };
  }

  /** RECIST-like call: burden is a volume, RECIST measures diameters, so convert with a cube root. */
  function recist(volRatio) {
    const diam = Math.cbrt(volRatio) - 1;
    if (volRatio < 1e-4) return { call: "CR", label: "Complete response", diam };
    if (diam <= -0.30) return { call: "PR", label: "Partial response", diam };
    if (diam >= 0.20) return { call: "PD", label: "Progressive disease", diam };
    return { call: "SD", label: "Stable disease", diam };
  }

  /** Simulated blood TCRβ-seq: multinomial draw of `reads` over infused clones + a diverse endogenous pool. */
  function bloodSample(tl, clones, day, reads, R) {
    const i = tl.t.indexOf(day);
    const sizes = tl.clones.map((v) => v[i]), infused = sizes.reduce((a, b) => a + b, 0), endo = tl.endo[i];
    const total = infused + endo;
    let detected = 0, top = 0, fromInfused = 0, remaining = reads, pLeft = 1;
    sizes.forEach((n) => { // sequential binomial draws = multinomial
      const p = n / total; const k = remaining > 0 && pLeft > 0 ? binom(R, remaining, Math.min(1, p / pLeft)) : 0;
      remaining -= k; pLeft -= p; fromInfused += k; if (k > 0) detected++; top = Math.max(top, k);
    });
    return { day, reads, persistence: fromInfused / reads, detected, topClone: top / reads, trueInfusedFrac: infused / total };
  }
  function binom(R, n, p) { // normal / poisson approximations are fine at these depths
    if (p <= 0) return 0; if (p >= 1) return n;
    const mu = n * p; if (mu < 20) return Math.min(n, poisson(R, mu));
    return Math.min(n, Math.max(0, Math.round(mu + Math.sqrt(mu * (1 - p)) * normal(R))));
  }

  function readouts(tl, clones, R) {
    const at = (d) => tl.t.indexOf(d);
    const wk12 = tl.tumor[at(84)] / tl.tumor[0], endRatio = tl.tumor[tl.tumor.length - 1] / tl.tumor[0];
    const nadir = Math.min(...tl.tumor), nadirDay = tl.t[tl.tumor.indexOf(nadir)];
    const end = tl.tumor[tl.tumor.length - 1];
    const relapse = nadir < 0.7 && end > Math.max(nadir * 10, 0.05);
    return {
      wk12: { ratio: wk12, ...recist(wk12) },
      best: { ratio: nadir, day: nadirDay, ...recist(nadir) },
      relapse,
      final: { day: tl.t[tl.t.length - 1], ratio: endRatio, ...recist(endRatio) },
      blood: [7, 30, 90, 180, 365].filter((d) => tl.t.includes(d)).map((d) => bloodSample(tl, clones, d, 20000, R)),
      survivors: tl.clones.filter((v) => v[v.length - 1] > 0).length,
    };
  }

  const api = { simulate, DEFAULTS, PARAMS, recist };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.SimModel = api;
})(typeof window !== "undefined" ? window : globalThis);
