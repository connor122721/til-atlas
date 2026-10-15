// Clone simulator UI. The model lives in sim-model.js (pure, testable in Node).
// Default input is REAL: tumor CD8 repertoires from Wu et al. 2020 (data/sim/repertoires.json).
(async function () {
  const { css, fmt, pct, esc } = TIL;
  const $ = (s) => document.querySelector(s);
  let seed = 7, last = null;
  const REPS = (await TIL.load("data/sim/repertoires.json")).data;

  const PRESETS = {
    "Responder-like": { stem: 0.35, reactive: 0.4, avidity: 1.2, clonality: 1.1 },
    "Stem-poor product": { stem: 0.02, reactive: 0.45, avidity: 1.2, clonality: 1.1 },
    "Mostly bystanders": { stem: 0.08, reactive: 0.05, avidity: 0.6, clonality: 1.1 },
  };
  const CLASS_CUSTOM = {
    high: { label: "Tumor-reactive, high avidity", color: "--s1" },
    low: { label: "Tumor-reactive, low avidity", color: "--s3" },
    bystander: { label: "Not tumor-reactive", color: "--s2" },
  };
  const CLASS_REAL = { high: { label: "Flagged tumor-reactive", color: "--s1" }, bystander: { label: "Not flagged", color: "--s2" } };
  const CLASS = () => (isReal() ? CLASS_REAL : CLASS_CUSTOM);
  const ICON = {
    CR: '<svg viewBox="0 0 20 20"><circle cx="10" cy="10" r="9" fill="var(--good)"/><path d="M6 10.5l2.5 2.5L14 7.5" stroke="#fff" stroke-width="2" fill="none"/></svg>',
    PR: '<svg viewBox="0 0 20 20"><circle cx="10" cy="10" r="9" fill="var(--good)"/><path d="M10 5v9M6.5 10.5 10 14l3.5-3.5" stroke="#fff" stroke-width="2" fill="none"/></svg>',
    SD: '<svg viewBox="0 0 20 20"><circle cx="10" cy="10" r="9" fill="var(--warning)"/><path d="M6 10h8" stroke="#000" stroke-width="2"/></svg>',
    PD: '<svg viewBox="0 0 20 20"><circle cx="10" cy="10" r="9" fill="var(--critical)"/><path d="M10 15V6M6.5 9.5 10 6l3.5 3.5" stroke="#fff" stroke-width="2" fill="none"/></svg>',
  };
  // Published TIL trial results for comparison (see page footnote for sources)
  const BENCH = [
    { name: "TIL arm, phase 3 (Rohaan 2022, NEJM)", n: 84, orr: 0.49, cr: 0.20 },
    { name: "Ipilimumab arm, same trial", n: 84, orr: 0.21, cr: 0.07 },
    { name: "Lifileucel, C-144-01 pooled (JITC 2022)", n: 153, orr: 48 / 153, cr: 8 / 153 },
  ];

  // ---------- source selector: real repertoires first, then a custom what-if product
  const reps = REPS ? REPS.patients : [];
  $("#source").innerHTML = reps.map((r, i) => `<option value="${i}">${esc(r.patient)} · ${esc(r.cancer)} tumor (Wu 2020)</option>`).join("") +
    `<option value="custom">Custom what-if product (sliders)</option>`;
  if (!reps.length) $("#source").value = "custom";
  else $("#source").value = String(Math.max(0, reps.findIndex((r) => r.patient === "Lung2")));
  const isReal = () => $("#source").value !== "custom";
  const rep = () => reps[+$("#source").value];
  const reactKey = () => $("#reactDef").value;

  const ctrl = {
    cells: (v) => `${(10 ** v / 1e9).toFixed(v >= 10 ? 0 : 1)}×10⁹`, stem: (v) => pct(v), reactive: (v) => pct(v),
    avidity: (v) => (+v).toFixed(2), clonality: (v) => (+v).toFixed(2), il2: (v) => (v == 0 ? "none" : `${v} d`),
    growth: (v) => `${v} d`, resist: (v) => `${(+v).toFixed(1)}×`,
  };
  const val = (id) => (id === "lymphodep" ? $("#" + id).checked : +$("#" + id).value);
  function syncOutputs() { Object.entries(ctrl).forEach(([id, f]) => ($(`#${id}-o`).textContent = f(val(id)))); }
  function realClones(r) { return r.clones.map((c) => ({ ...c, reactive: !!c[reactKey()] })); }
  function settings() {
    const base = { productCells: 10 ** val("cells") / 1e6, avidity: val("avidity"), lymphodep: val("lymphodep"), il2: val("il2") > 0,
      il2Days: val("il2"), tumorGrowth: Math.LN2 / val("growth"), resistance: val("resist") };
    if (isReal()) return { ...base, repertoire: realClones(rep()) };
    return { ...base, stemFrac: val("stem"), reactiveFrac: val("reactive"), clonality: val("clonality") };
  }
  function stemOfReactive(clones) {
    const rx = clones.filter((c) => c.reactive), tot = d3.sum(rx, (c) => c.S + c.E + c.X);
    return tot ? d3.sum(rx, (c) => c.S) / tot : null;
  }
  function describeProduct() {
    $("#realOpts").hidden = !isReal(); $("#customOpts").hidden = isReal();
    $("#cohort").textContent = isReal() ? `Infuse all ${reps.length} real repertoires` : "Run 40 hypothetical products";
    if (!isReal()) return;
    const r = rep(), cl = realClones(r), tot = d3.sum(cl, (c) => c.S + c.E + c.X), rx = cl.filter((c) => c.reactive);
    const sr = stemOfReactive(cl);
    $("#productNote").innerHTML = `${r.clones.length} largest of ${fmt(r.clonesTotal)} tumor CD8⁺ clones (${fmt(tot)} of ${fmt(r.tumorCD8)} cells).
      Stem-like ${pct(d3.sum(cl, (c) => c.S) / tot)}, terminal ${pct(d3.sum(cl, (c) => c.X) / tot)} of cells.
      <b>${rx.length}</b> clones flagged tumor-reactive${rx.length ? `, ${pct(sr)} of their cells stem-like` : " — none met the rule, so no tumor-reactive cells are infused"}.`;
  }

  $("#presets").innerHTML = Object.keys(PRESETS).map((k) => `<button type="button" class="chip" data-p="${esc(k)}">${esc(k)}</button>`).join("");
  $("#presets").addEventListener("click", (e) => {
    const b = e.target.closest("[data-p]"); if (!b) return;
    Object.entries(PRESETS[b.dataset.p]).forEach(([k, v]) => ($("#" + k).value = v));
    document.querySelectorAll("#presets .chip").forEach((c) => c.setAttribute("aria-pressed", c === b));
    run();
  });
  let timer;
  document.querySelectorAll(".controls input").forEach((el) => el.addEventListener("input", () => {
    document.querySelectorAll("#presets .chip").forEach((c) => c.setAttribute("aria-pressed", "false"));
    syncOutputs(); clearTimeout(timer); timer = setTimeout(run, 120);
  }));
  ["#source", "#reactDef"].forEach((s) => $(s).addEventListener("change", () => { cohort = null; $("#cohortCard").hidden = true; run(); }));
  $("#reseed").addEventListener("click", () => { seed = Math.floor(Math.random() * 1e6); run(); });

  // ---------- single run
  function run() {
    syncOutputs(); describeProduct();
    last = SimModel.simulate({ ...settings(), seed });
    $("#seedNote").textContent = `Random seed ${seed}. ${isReal() ? "Clone sizes and states are measured; divisions, deaths and avidities are drawn at random." : "Every value in a custom product is hypothetical."}`;
    drawTiles(last.readouts); drawMuller(); drawStates(); drawTumor(); drawBlood(last.readouts);
  }

  function drawTiles(r) {
    const d30 = r.blood.find((b) => b.day === 30);
    const tile = (label, value, note) => `<div class="card stat"><div class="label">${label}</div><div class="value">${value}</div><div class="note">${note}</div></div>`;
    const resp = (x) => `<span class="outcome">${ICON[x.call]}${x.call}</span>`;
    $("#tiles").innerHTML = [
      tile("Best response", resp(r.best), `${r.best.label}, day ${r.best.day} (${diamTxt(r.best.diam)} diameter)`),
      tile("Week 12", resp(r.wk12), `${r.wk12.label} (${diamTxt(r.wk12.diam)})`),
      tile("Day 365", resp(r.final), r.relapse ? "Responded, then relapsed" : r.final.label),
      tile("Persistence, day 30", pct(d30.persistence), `top clone ${pct(d30.topClone, 1)} of blood TCRs${d30.topClone > 0.05 ? " — immunodominant" : ""}`),
    ].join("");
  }
  const diamTxt = (d) => (d <= -0.999 ? "−100%" : `${d > 0 ? "+" : ""}${Math.round(d * 100)}%`);

  // ---------- Muller-style plot of clone shares
  function drawMuller() {
    const el = $("#muller"); el.innerHTML = "";
    const { timeline: tl, clones } = last, CL = CLASS();
    const idx = tl.t.map((_, i) => i).filter((i) => i % 2 === 0);
    const order = clones.map((c, i) => i).sort((a, b) => {
      const rank = { high: 0, low: 1, bystander: 2 };
      return rank[clones[a].cls] - rank[clones[b].cls] || clones[b].infused - clones[a].infused;
    });
    const rows = idx.map((i) => { const o = { t: tl.t[i] }; order.forEach((k) => (o[k] = tl.clones[k][i])); return o; });
    const stack = d3.stack().keys(order).offset(d3.stackOffsetExpand)(rows);
    const W = el.clientWidth || 700, H = 300, m = { t: 6, r: 8, b: 26, l: 40 };
    const x = d3.scaleLinear().domain([0, tl.t[tl.t.length - 1]]).range([m.l, W - m.r]);
    const y = d3.scaleLinear().domain([0, 1]).range([H - m.b, m.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img")
      .attr("aria-label", "Stacked area chart of each infused clone's share over time");
    const area = d3.area().x((d) => x(d.data.t)).y0((d) => y(isNaN(d[0]) ? 0 : d[0])).y1((d) => y(isNaN(d[1]) ? 0 : d[1]));
    const maxInf = d3.max(clones, (c) => c.infused);
    svg.append("g").selectAll("path").data(stack).join("path")
      .attr("d", area).attr("fill", (s) => css(CL[clones[s.key].cls].color))
      .attr("fill-opacity", (s) => 0.55 + 0.45 * Math.sqrt(clones[s.key].infused / maxInf))
      .attr("stroke", css("--surface-1")).attr("stroke-width", 0.6);
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`)
      .call(d3.axisBottom(x).ticks(W > 500 ? 8 : 4).tickFormat((d) => `d${d}`).tickSize(0).tickPadding(7)).call((g) => g.select(".domain").remove());
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`)
      .call(d3.axisLeft(y).ticks(4, "%").tickSize(0).tickPadding(6)).call((g) => g.select(".domain").remove());
    const cross = svg.append("line").attr("y1", m.t).attr("y2", H - m.b).attr("stroke", css("--text-muted")).attr("stroke-width", 1).style("opacity", 0);
    svg.append("rect").attr("x", m.l).attr("y", m.t).attr("width", W - m.l - m.r).attr("height", H - m.t - m.b).attr("fill", "transparent")
      .on("mousemove", (ev) => {
        const [px, py] = d3.pointer(ev);
        const j = Math.max(0, Math.min(rows.length - 1, Math.round((x.invert(px) / tl.t[tl.t.length - 1]) * (rows.length - 1))));
        const fy = y.invert(py), layer = stack.find((s) => s[j][0] <= fy && fy <= s[j][1]);
        cross.attr("x1", x(rows[j].t)).attr("x2", x(rows[j].t)).style("opacity", 1);
        if (!layer) return TIL.hideTip();
        const c = clones[layer.key], n = rows[j][layer.key];
        TIL.showTip(`<div class="tt-title">${c.label ? `<span class="mono">${esc(c.label)}</span>` : `Clone #${c.id}`} · day ${rows[j].t}</div>
          <div class="row"><i style="background:${css(CL[c.cls].color)}"></i>${CL[c.cls].label}</div>
          <div class="row">Share of infused cells<b>${pct(layer[j][1] - layer[j][0], 1)}</b></div>
          <div class="row">Cells (×10⁶)<b>${fmt(n)}</b></div>
          <div class="row">At infusion<b>${fmt(c.infused)} · ${pct(c.stem0)} stem-like</b></div>
          ${c.label ? `<div class="row">Measured in tumor<b>${pct(c.stem0)} stem-like${c.inNAT ? ` · ${c.inNAT} cells in NAT` : ""}${c.inBlood ? " · in blood" : ""}</b></div>` : ""}
          ${c.reactive ? `<div class="row">Avidity (drawn)<b>${c.avidity.toFixed(2)}</b></div>` : ""}`, ev);
      })
      .on("mouseleave", () => { TIL.hideTip(); cross.style("opacity", 0); });
    $("#mullerLegend").innerHTML = Object.values(CL).map((c) => `<span><i style="background:${css(c.color)}"></i>${c.label}</span>`).join("");
  }

  // ---------- log-scale line chart with crosshair tooltip
  function lineChart(el, series, { yMin, refs = [], label, tickFmt = d3.format("~s") }) {
    el.innerHTML = "";
    const t = last.timeline.t, W = el.clientWidth || 400, H = 230, m = { t: 8, r: 64, b: 26, l: 44 };
    const x = d3.scaleLinear().domain([0, t[t.length - 1]]).range([m.l, W - m.r]);
    const yMax = d3.max(series, (s) => d3.max(s.values));
    const y = d3.scaleLog().domain([yMin, Math.max(yMax * 1.5, yMin * 10)]).range([H - m.b, m.t]).clamp(true);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", label);
    const ticks = y.ticks().filter((v) => Math.log10(v) % 1 === 0);
    svg.append("g").attr("class", "gridline").attr("transform", `translate(${m.l},0)`)
      .call(d3.axisLeft(y).tickValues(ticks).tickSize(-(W - m.l - m.r)).tickFormat(""));
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`)
      .call(d3.axisLeft(y).tickValues(ticks).tickFormat(tickFmt).tickSize(0).tickPadding(6)).call((g) => g.select(".domain").remove());
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`)
      .call(d3.axisBottom(x).ticks(W > 420 ? 6 : 3).tickFormat((d) => `d${d}`).tickSize(0).tickPadding(7)).call((g) => g.select(".domain").remove());
    refs.forEach((r) => {
      svg.append("line").attr("x1", m.l).attr("x2", W - m.r).attr("y1", y(r.y)).attr("y2", y(r.y))
        .attr("stroke", css("--text-muted")).attr("stroke-dasharray", "4 4");
      svg.append("text").attr("x", W - m.r + 6).attr("y", y(r.y)).attr("dy", "0.35em").attr("font-size", 11).attr("fill", css("--text-muted")).text(r.label);
    });
    const line = d3.line().defined((v) => v > 0).x((_, i) => x(t[i])).y((v) => y(Math.max(v, yMin)));
    series.forEach((s) => {
      svg.append("path").attr("d", line(s.values)).attr("fill", "none").attr("stroke", css(s.color)).attr("stroke-width", 2)
        .attr("stroke-dasharray", s.dash || null);
      const lastV = s.values[s.values.length - 1];
      if (series.length > 1 && lastV > yMin)
        svg.append("text").attr("x", W - m.r + 6).attr("y", y(lastV)).attr("dy", "0.35em").attr("font-size", 11.5)
          .attr("fill", css("--text-secondary")).text(s.short || s.name);
    });
    const cross = svg.append("line").attr("y1", m.t).attr("y2", H - m.b).attr("stroke", css("--text-muted")).style("opacity", 0);
    const dots = svg.selectAll("circle.hov").data(series).join("circle").attr("class", "hov").attr("r", 4)
      .attr("fill", (s) => css(s.color)).attr("stroke", css("--surface-1")).attr("stroke-width", 2).style("opacity", 0);
    svg.append("rect").attr("x", m.l).attr("y", m.t).attr("width", W - m.l - m.r).attr("height", H - m.t - m.b).attr("fill", "transparent")
      .on("mousemove", (ev) => {
        const i = Math.max(0, Math.min(t.length - 1, Math.round(x.invert(d3.pointer(ev)[0]) - t[0])));
        cross.attr("x1", x(t[i])).attr("x2", x(t[i])).style("opacity", 1);
        dots.attr("cx", x(t[i])).attr("cy", (s) => y(Math.max(s.values[i], yMin))).style("opacity", (s) => (s.values[i] > 0 ? 1 : 0));
        TIL.showTip(`<div class="tt-title">Day ${t[i]}</div>` + series.map((s) =>
          `<div class="row"><i style="background:${css(s.color)}"></i>${s.name}<b>${s.fmt ? s.fmt(s.values[i]) : fmt(s.values[i])}</b></div>`).join(""), ev);
      })
      .on("mouseleave", () => { TIL.hideTip(); cross.style("opacity", 0); dots.style("opacity", 0); });
  }
  function drawStates() {
    const tl = last.timeline;
    const series = [
      { name: "Stem-like (S)", short: "S", values: tl.S, color: "--s1" },
      { name: "Effector (E)", short: "E", values: tl.E, color: "--s2" },
      { name: "Exhausted (X)", short: "X", values: tl.X, color: "--s3" },
      { name: "Endogenous T cells", short: "Endog.", values: tl.endo, color: "--text-muted", dash: "5 4" },
    ];
    $("#stateLegend").innerHTML = series.map((s) => `<span><i style="background:${css(s.color)}"></i>${s.name}</span>`).join("");
    lineChart($("#states"), series, { yMin: 1, label: "Infused T cells by phenotype state over time, log scale" });
  }
  function drawTumor() {
    const v = last.timeline.tumor;
    lineChart($("#tumor"), [{ name: "Tumor burden", values: v, color: "--s8", fmt: (x) => (x === 0 ? "cleared" : `${(x * 100).toPrecision(2)}% of baseline`) }],
      { yMin: 1e-4, label: "Tumor burden relative to baseline, log scale", tickFmt: (v) => `${+(v * 100).toPrecision(1)}%`,
        refs: [{ y: 0.343, label: "PR −30%" }, { y: 1.728, label: "PD +20%" }] });
  }
  function drawBlood(r) {
    $("#blood tbody").innerHTML = r.blood.map((b) => `<tr><td>Day ${b.day}</td><td class="num">${pct(b.persistence, 1)}</td>
      <td class="num">${b.detected} / ${last.clones.length}</td><td class="num">${pct(b.topClone, 1)}</td><td class="num">${pct(b.trueInfusedFrac, 1)}</td></tr>`).join("");
  }

  // ---------- virtual trial
  let cohort = null;
  const SEEDS = 5;
  $("#cohort").addEventListener("click", async () => {
    const btn = $("#cohort"); btn.disabled = true;
    const base = settings(), out = [];
    if (isReal()) {
      for (let i = 0; i < reps.length; i++) {
        const cl = realClones(reps[i]), runs = [];
        for (let k = 0; k < SEEDS; k++) runs.push(SimModel.simulate({ ...base, repertoire: cl, seed: 1000 + 37 * i + k, days: 180 }).readouts);
        const diams = runs.map((r) => r.best.diam).sort(d3.ascending), med = runs[runs.map((r) => r.best.diam).indexOf(d3.quantileSorted(diams, 0.5))] || runs[0];
        out.push({ i: i + 1, label: reps[i].patient, cancer: reps[i].cancer, runs, best: med.best, final: med.final, relapse: med.relapse,
          lo: diams[0], hi: diams[diams.length - 1], x: stemOfReactive(cl), nReactive: cl.filter((c) => c.reactive).length });
        btn.textContent = `Infusing ${reps[i].patient}…`; await new Promise((res) => setTimeout(res));
      }
    } else {
      for (let i = 0; i < 40; i++) {
        const rn = d3.randomNormal.source(d3.randomLcg(((i + 1) * 0.6180339887) % 1))(0, 1), jit = (sd) => Math.exp(sd * rn());
        const p = { stemFrac: Math.min(0.9, base.stemFrac * jit(0.8)), reactiveFrac: Math.min(0.9, base.reactiveFrac * jit(0.6)) };
        const r = SimModel.simulate({ ...base, ...p, seed: 1000 + i, days: 180, avidity: base.avidity * jit(0.7),
          tumorGrowth: base.tumorGrowth * jit(0.3), resistance: base.resistance * jit(0.8) }).readouts;
        out.push({ i: i + 1, label: `Hypothetical ${i + 1}`, runs: [r], best: r.best, final: r.final, relapse: r.relapse, x: p.stemFrac });
        if (i % 5 === 4) { btn.textContent = `Product ${i + 1}/40…`; await new Promise((res) => setTimeout(res)); }
      }
    }
    cohort = { real: isReal(), rows: out }; btn.disabled = false; describeProduct();
    $("#cohortCard").hidden = false; drawWaterfall(); $("#cohortCard").scrollIntoView({ behavior: "smooth", block: "nearest" });
  });
  function drawWaterfall() {
    if (!cohort) return;
    const el = $("#waterfall"); el.innerHTML = "";
    const data = cohort.rows.map((d) => ({ ...d, v: Math.max(-1, Math.min(1, d.best.diam)) })).sort((a, b) => b.v - a.v);
    const allRuns = cohort.rows.flatMap((d) => d.runs), n = (c) => allRuns.filter((r) => r.best.call === c).length;
    const orr = (n("CR") + n("PR")) / allRuns.length, cr = n("CR") / allRuns.length, rel = allRuns.filter((r) => r.relapse).length;
    $("#cohortTitle").textContent = cohort.real ? "Virtual trial: each real repertoire infused" : "Virtual trial: 40 hypothetical products";
    $("#cohortSub").textContent = cohort.real
      ? `${cohort.rows.length} real tumor repertoires × ${SEEDS} random draws = ${allRuns.length} virtual infusions. Objective response ${pct(orr)}, complete ${pct(cr)}, ` +
        `${rel} relapsed by day 180. Bars show each repertoire's median best change; whiskers span its ${SEEDS} draws.`
      : `Objective response ${pct(orr)}, complete ${pct(cr)}, ${rel} relapsed by day 180. Every product here is hypothetical.`;
    const cats = [["Response (CR/PR)", "--good"], ["Stable disease", "--warning"], ["Progression", "--critical"]];
    $("#cohortLegend").innerHTML = cats.map(([l, c]) => `<span><i style="background:${css(c)}"></i>${l}</span>`).join("");
    const W = el.clientWidth || 700, H = 250, m = { t: 10, r: 8, b: cohort.real ? 46 : 14, l: 44 };
    const x = d3.scaleBand().domain(data.map((d) => d.i)).range([m.l, W - m.r]).paddingInner(0.25);
    const y = d3.scaleLinear().domain([-1, 1]).range([H - m.b, m.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Waterfall plot of best response");
    svg.append("g").attr("class", "gridline").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([-1, -0.3, 0.2, 1]).tickSize(-(W - m.l - m.r)).tickFormat(""));
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`)
      .call(d3.axisLeft(y).tickValues([-1, -0.3, 0, 0.2, 1]).tickFormat((d) => `${d > 0 ? "+" : ""}${Math.round(d * 100)}%`).tickSize(0).tickPadding(6))
      .call((g) => g.select(".domain").remove());
    const bw = x.bandwidth(), r = Math.min(4, bw / 2);
    const path = (d) => { const x0 = x(d.i), y0 = y(0), y1 = y(d.v);
      if (Math.abs(y1 - y0) < 3) return `M${x0},${y0 - 1.5}h${bw}v3h${-bw}Z`;   // no change: keep a visible tick at 0
      return d.v > 0 ? TILCharts.barPath(x0, y1, bw, y0 - y1, r)
        : `M${x0},${y0}V${y1 - r}Q${x0},${y1} ${x0 + r},${y1}H${x0 + bw - r}Q${x0 + bw},${y1} ${x0 + bw},${y1 - r}V${y0}Z`; };
    const color = (d) => css(d.best.call === "CR" || d.best.call === "PR" ? "--good" : d.best.call === "SD" ? "--warning" : "--critical");
    svg.append("g").selectAll("path").data(data).join("path").attr("d", path).attr("fill", color);
    if (cohort.real) {
      svg.append("g").selectAll("line").data(data).join("line").attr("x1", (d) => x(d.i) + bw / 2).attr("x2", (d) => x(d.i) + bw / 2)
        .attr("y1", (d) => y(Math.max(-1, Math.min(1, d.lo)))).attr("y2", (d) => y(Math.max(-1, Math.min(1, d.hi)))).attr("stroke", css("--text-secondary")).attr("stroke-width", 1.5);
      svg.append("g").selectAll("text").data(data).join("text").attr("transform", (d) => `translate(${x(d.i) + bw / 2},${H - m.b + 8}) rotate(-40)`)
        .attr("text-anchor", "end").attr("font-size", 11).attr("fill", css("--text-secondary")).text((d) => d.label);
    }
    svg.append("line").attr("x1", m.l).attr("x2", W - m.r).attr("y1", y(0)).attr("y2", y(0)).attr("stroke", css("--axis"));
    svg.append("g").selectAll("rect").data(data).join("rect").attr("x", (d) => x(d.i)).attr("width", bw).attr("y", m.t).attr("height", H - m.t - m.b)
      .attr("fill", "transparent")
      .on("mousemove", (ev, d) => TIL.showTip(`<div class="tt-title">${esc(d.label)}${d.cancer ? ` · ${esc(d.cancer)}` : ""}</div>
        <div class="row">Best response (median draw)<b>${d.best.call} (${diamTxt(d.best.diam)})</b></div>
        <div class="row">Day 180<b>${d.final.call}${d.relapse ? ", relapsed" : ""}</b></div>
        ${d.runs.length > 1 ? `<div class="row">Responses in ${d.runs.length} draws<b>${d.runs.filter((r) => ["CR", "PR"].includes(r.best.call)).length}</b></div>` : ""}
        ${cohort.real && !d.nReactive ? `<div class="small muted">No clone met the tumor-reactive rule, so nothing in the product attacks the tumor.</div>` : ""}
        ${cohort.real ? `<div class="row">Flagged tumor-reactive clones<b>${d.nReactive}</b></div><div class="row">Stem-like share of their cells<b>${d.x == null ? "—" : pct(d.x)}</b></div>`
          : `<div class="row">Product stem-like share<b>${pct(d.x)}</b></div>`}`, ev))
      .on("mouseleave", TIL.hideTip);
    drawBench(orr, cr, allRuns.length);
    drawScatter(data.filter((d) => d.x != null), color);
  }

  function drawBench(orr, cr, n) {
    const el = $("#bench"); el.innerHTML = "";
    const rows = [{ name: cohort.real ? "This virtual trial (real repertoires)" : "This virtual trial (hypothetical products)", n, orr, cr, virtual: true }, ...BENCH];
    const W = el.clientWidth || 700, rowH = 30, lab = Math.min(270, W * 0.45), H = rows.length * rowH + 26;
    const x = d3.scaleLinear().domain([0, 1]).range([lab, W - 50]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Response rates compared with published TIL trials");
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - 22})`).call(d3.axisBottom(x).ticks(5, "%").tickSize(0).tickPadding(6)).call((g) => g.select(".domain").remove());
    rows.forEach((r, i) => {
      const y0 = i * rowH;
      svg.append("text").attr("x", lab - 8).attr("y", y0 + rowH / 2).attr("dy", "0.35em").attr("text-anchor", "end").attr("font-size", 12)
        .attr("font-weight", r.virtual ? 700 : 400).attr("fill", css(r.virtual ? "--text-primary" : "--text-secondary")).text(r.name);
      svg.append("rect").attr("x", x(0)).attr("y", y0 + 6).attr("width", x(r.cr) - x(0)).attr("height", rowH - 12).attr("fill", css("--s1"));
      svg.append("path").attr("d", TILCharts.barPath(x(r.cr), y0 + 6, Math.max(0, x(r.orr) - x(r.cr)), rowH - 12, 4, true)).attr("fill", css("--q300"));
      svg.append("text").attr("x", x(r.orr) + 6).attr("y", y0 + rowH / 2).attr("dy", "0.35em").attr("font-size", 11.5).attr("fill", css("--text-secondary"))
        .style("font-variant-numeric", "tabular-nums").text(`${pct(r.orr)}`);
      svg.append("rect").attr("x", 0).attr("y", y0).attr("width", W).attr("height", rowH).attr("fill", "transparent")
        .on("mousemove", (ev) => TIL.showTip(`<div class="tt-title">${esc(r.name)}</div><div class="row">Objective response<b>${pct(r.orr, 1)}</b></div>
          <div class="row">Complete response<b>${pct(r.cr, 1)}</b></div><div class="row">${r.virtual ? "Virtual infusions" : "Patients"}<b>${r.n}</b></div>`, ev))
        .on("mouseleave", TIL.hideTip);
    });
    el.insertAdjacentHTML("afterbegin", `<div class="legend"><span><i style="background:${css("--s1")}"></i>Complete response</span><span><i style="background:${css("--q300")}"></i>Partial response</span></div>`);
  }

  // Krishna inside the simulator: does product stemness predict the best response?
  function drawScatter(data, color) {
    const el = $("#scatter"); el.innerHTML = "";
    const W = el.clientWidth || 600, H = 240, m = { t: 10, r: 12, b: 36, l: 44 };
    if (!data.length) { el.innerHTML = `<p class="muted small">No repertoire had flagged tumor-reactive clones.</p>`; return; }
    const x = d3.scaleLinear().domain([0, Math.max(0.05, d3.max(data, (d) => d.x) * 1.05)]).nice().range([m.l, W - m.r]);
    const y = d3.scaleLinear().domain([-1, 1]).range([H - m.b, m.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Scatter of product stem-like share against best response");
    svg.append("g").attr("class", "gridline").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([-1, -0.3, 0.2, 1]).tickSize(-(W - m.l - m.r)).tickFormat(""));
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`)
      .call(d3.axisLeft(y).tickValues([-1, -0.3, 0, 0.2, 1]).tickFormat((d) => `${d > 0 ? "+" : ""}${Math.round(d * 100)}%`).tickSize(0).tickPadding(6)).call((g) => g.select(".domain").remove());
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`)
      .call(d3.axisBottom(x).ticks(6, "%").tickSize(0).tickPadding(7)).call((g) => g.select(".domain").remove());
    svg.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 4).attr("text-anchor", "middle").attr("font-size", 12)
      .attr("fill", css("--text-secondary")).text(cohort.real ? "Stem-like share of the flagged tumor-reactive clones, measured in the tumor" : "Stem-like share of the hypothetical product");
    svg.append("g").selectAll("circle").data(data).join("circle").attr("cx", (d) => x(d.x)).attr("cy", (d) => y(d.v)).attr("r", 5.5)
      .attr("fill", color).attr("stroke", css("--surface-1")).attr("stroke-width", 2)
      .on("mousemove", (ev, d) => TIL.showTip(`<div class="tt-title">${esc(d.label)}</div>
        <div class="row">Stem-like share<b>${pct(d.x)}</b></div><div class="row">Best response<b>${d.best.call} (${diamTxt(d.best.diam)})</b></div>
        <div class="row">Day 180<b>${d.final.call}${d.relapse ? ", relapsed" : ""}</b></div>`, ev))
      .on("mouseleave", TIL.hideTip);
    const rho = spearman(data.map((d) => d.x), data.map((d) => d.v));
    const medx = d3.median(data, (e) => e.x);
    $("#scatterNote").textContent = `Spearman ρ = ${rho.toFixed(2)} between stemness and best tumor change (negative = more stem-like, more shrinkage), n = ${data.length}. ` +
      `Relapse by day 180: ${pct(d3.mean(data.filter((d) => d.x < medx), (d) => +d.relapse))} below the median stemness vs ` +
      `${pct(d3.mean(data.filter((d) => d.x >= medx), (d) => +d.relapse))} at or above. ` +
      (cohort.real ? "The x-axis is measured; the response is the model's." : "");
  }
  function spearman(a, b) {
    const rank = (v) => { const o = v.map((x, i) => [x, i]).sort((p, q) => p[0] - q[0]), r = Array(v.length); o.forEach(([, i], k) => (r[i] = k)); return r; };
    const ra = rank(a), rb = rank(b), ma = d3.mean(ra), mb = d3.mean(rb);
    return d3.sum(ra, (x, i) => (x - ma) * (rb[i] - mb)) / Math.sqrt(d3.sum(ra, (x) => (x - ma) ** 2) * d3.sum(rb, (x) => (x - mb) ** 2));
  }

  let rt; window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(() => { if (last) { drawMuller(); drawStates(); drawTumor(); drawWaterfall(); } }, 150); });
  window.addEventListener("themechange", () => { if (last) { drawMuller(); drawStates(); drawTumor(); drawWaterfall(); } });
  $("#reactDef").value = "enriched";
  run();
})();
