// Clone simulator UI. The model itself lives in sim-model.js (pure, testable in Node).
(function () {
  const { css, fmt, pct, esc } = TIL;
  const $ = (s) => document.querySelector(s);
  let seed = 7, last = null;

  const PRESETS = {
    "Responder-like": { cells: 10.5, stem: 0.35, reactive: 0.4, avidity: 1.2, clonality: 1.1, lymphodep: true, il2: 3, growth: 20, resist: 1 },
    "Stem-poor product": { cells: 10.5, stem: 0.02, reactive: 0.45, avidity: 1.2, clonality: 1.1, lymphodep: true, il2: 3, growth: 20, resist: 1 },
    "Blood-borne bystanders": { cells: 10.5, stem: 0.08, reactive: 0.05, avidity: 0.6, clonality: 1.1, lymphodep: true, il2: 3, growth: 20, resist: 1 },
    "No lymphodepletion": { cells: 10.5, stem: 0.35, reactive: 0.4, avidity: 1.2, clonality: 1.1, lymphodep: false, il2: 3, growth: 20, resist: 1 },
    "Fast, resistant tumor": { cells: 10.5, stem: 0.3, reactive: 0.35, avidity: 1, clonality: 1.1, lymphodep: true, il2: 3, growth: 10, resist: 2.5 },
  };
  const CLASS = {
    high: { label: "Tumor-reactive, high avidity", color: "--s1" },
    low: { label: "Tumor-reactive, low avidity", color: "--s3" },
    bystander: { label: "Bystander (not tumor-reactive)", color: "--s2" },
  };
  const OUT = { CR: "--good", PR: "--good", SD: "--warning", PD: "--critical" };
  const ICON = {
    CR: '<svg viewBox="0 0 20 20"><circle cx="10" cy="10" r="9" fill="var(--good)"/><path d="M6 10.5l2.5 2.5L14 7.5" stroke="#fff" stroke-width="2" fill="none"/></svg>',
    PR: '<svg viewBox="0 0 20 20"><circle cx="10" cy="10" r="9" fill="var(--good)"/><path d="M10 5v9M6.5 10.5 10 14l3.5-3.5" stroke="#fff" stroke-width="2" fill="none"/></svg>',
    SD: '<svg viewBox="0 0 20 20"><circle cx="10" cy="10" r="9" fill="var(--warning)"/><path d="M6 10h8" stroke="#000" stroke-width="2"/></svg>',
    PD: '<svg viewBox="0 0 20 20"><circle cx="10" cy="10" r="9" fill="var(--critical)"/><path d="M10 15V6M6.5 9.5 10 6l3.5 3.5" stroke="#fff" stroke-width="2" fill="none"/></svg>',
  };

  // ---------- controls
  const ctrl = {
    cells: (v) => `${(10 ** v / 1e9).toFixed(v >= 10 ? 0 : 1)}×10⁹`, stem: (v) => pct(v), reactive: (v) => pct(v),
    avidity: (v) => (+v).toFixed(2), clonality: (v) => (+v).toFixed(2), il2: (v) => (v == 0 ? "none" : `${v} d`),
    growth: (v) => `${v} d`, resist: (v) => `${(+v).toFixed(1)}×`,
  };
  const val = (id) => (id === "lymphodep" ? $("#" + id).checked : +$("#" + id).value);
  function syncOutputs() { Object.entries(ctrl).forEach(([id, f]) => ($(`#${id}-o`).textContent = f(val(id)))); }
  function settings() {
    return { productCells: 10 ** val("cells") / 1e6, stemFrac: val("stem"), reactiveFrac: val("reactive"), avidity: val("avidity"),
      clonality: val("clonality"), lymphodep: val("lymphodep"), il2: val("il2") > 0, il2Days: val("il2"),
      tumorGrowth: Math.LN2 / val("growth"), resistance: val("resist") };
  }
  $("#presets").innerHTML = Object.keys(PRESETS).map((k) => `<button type="button" class="chip" data-p="${esc(k)}">${esc(k)}</button>`).join("");
  $("#presets").addEventListener("click", (e) => {
    const b = e.target.closest("[data-p]"); if (!b) return;
    const p = PRESETS[b.dataset.p];
    Object.entries(p).forEach(([k, v]) => (k === "lymphodep" ? ($("#lymphodep").checked = v) : ($("#" + k).value = v)));
    document.querySelectorAll("#presets .chip").forEach((c) => c.setAttribute("aria-pressed", c === b));
    run();
  });
  let timer;
  document.querySelectorAll(".controls input").forEach((el) => el.addEventListener("input", () => {
    document.querySelectorAll("#presets .chip").forEach((c) => c.setAttribute("aria-pressed", "false"));
    syncOutputs(); clearTimeout(timer); timer = setTimeout(run, 120);
  }));
  $("#reseed").addEventListener("click", () => { seed = Math.floor(Math.random() * 1e6); run(); });

  // ---------- single run
  function run() {
    syncOutputs();
    last = SimModel.simulate({ ...settings(), seed });
    $("#seedNote").textContent = `Random seed ${seed} — same seed + same settings = same patient.`;
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
    const { timeline: tl, clones } = last;
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
      .attr("d", area).attr("fill", (s) => css(CLASS[clones[s.key].cls].color))
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
        TIL.showTip(`<div class="tt-title">Clone #${c.id} · day ${rows[j].t}</div>
          <div class="row"><i style="background:${css(CLASS[c.cls].color)}"></i>${CLASS[c.cls].label}</div>
          <div class="row">Share of infused cells<b>${pct(layer[j][1] - layer[j][0], 1)}</b></div>
          <div class="row">Cells (×10⁶)<b>${fmt(n)}</b></div>
          <div class="row">At infusion<b>${fmt(c.infused)} · ${pct(c.stem0)} stem-like</b></div>
          ${c.reactive ? `<div class="row">Avidity<b>${c.avidity.toFixed(2)}</b></div>` : ""}`, ev);
      })
      .on("mouseleave", () => { TIL.hideTip(); cross.style("opacity", 0); });
    $("#mullerLegend").innerHTML = Object.values(CLASS).map((c) => `<span><i style="background:${css(c.color)}"></i>${c.label}</span>`).join("");
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
  $("#cohort").addEventListener("click", async () => {
    const btn = $("#cohort"); btn.disabled = true;
    const base = settings(), out = [];
    for (let i = 0; i < 40; i++) {
      // patients differ in their tumors AND their products: each manufactured product has its own
      // stem-like share and reactive fraction (log-normal around the slider values), as real products do
      const rn = d3.randomNormal.source(d3.randomLcg(((i + 1) * 0.6180339887) % 1))(0, 1), jit = (sd) => Math.exp(sd * rn());
      const p = { stemFrac: Math.min(0.9, base.stemFrac * jit(0.8)), reactiveFrac: Math.min(0.9, base.reactiveFrac * jit(0.6)) };
      const r = SimModel.simulate({ ...base, ...p, seed: 1000 + i, days: 180,
        avidity: base.avidity * jit(0.7), tumorGrowth: base.tumorGrowth * jit(0.3), resistance: base.resistance * jit(0.8) }).readouts;
      out.push({ i: i + 1, best: r.best, relapse: r.relapse, final: r.final, ...p });
      if (i % 5 === 4) { btn.textContent = `Treating patient ${i + 1}/40…`; await new Promise((res) => setTimeout(res)); }
    }
    cohort = out; btn.disabled = false; btn.textContent = "Run a 40-patient trial";
    $("#cohortCard").hidden = false; drawWaterfall(); $("#cohortCard").scrollIntoView({ behavior: "smooth", block: "nearest" });
  });
  function drawWaterfall() {
    if (!cohort) return;
    const el = $("#waterfall"); el.innerHTML = "";
    const data = [...cohort].map((d) => ({ ...d, v: Math.max(-1, Math.min(1, d.best.diam)) })).sort((a, b) => b.v - a.v);
    const n = (c) => data.filter((d) => d.best.call === c).length;
    const orr = n("CR") + n("PR");
    $("#cohortSub").textContent = `Objective response ${orr}/40 (${pct(orr / 40)}), complete ${n("CR")}/40, ` +
      `${data.filter((d) => d.relapse).length} relapsed by day 180. Best change in tumor diameter, one bar per virtual patient.`;
    const cats = [["Response (CR/PR)", "--good"], ["Stable disease", "--warning"], ["Progression", "--critical"]];
    $("#cohortLegend").innerHTML = cats.map(([l, c]) => `<span><i style="background:${css(c)}"></i>${l}</span>`).join("");
    const W = el.clientWidth || 700, H = 240, m = { t: 10, r: 8, b: 14, l: 44 };
    const x = d3.scaleBand().domain(data.map((d) => d.i)).range([m.l, W - m.r]).paddingInner(0.2);
    const y = d3.scaleLinear().domain([-1, 1]).range([H - m.b, m.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Waterfall plot of best response per virtual patient");
    svg.append("g").attr("class", "gridline").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([-1, -0.3, 0.2, 1]).tickSize(-(W - m.l - m.r)).tickFormat(""));
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`)
      .call(d3.axisLeft(y).tickValues([-1, -0.3, 0, 0.2, 1]).tickFormat((d) => `${d > 0 ? "+" : ""}${Math.round(d * 100)}%`).tickSize(0).tickPadding(6))
      .call((g) => g.select(".domain").remove());
    const bw = x.bandwidth(), r = Math.min(4, bw / 2);
    const path = (d) => { const x0 = x(d.i), y0 = y(0), y1 = y(d.v); if (Math.abs(y1 - y0) < 1) return "";
      return d.v > 0 ? TILCharts.barPath(x0, y1, bw, y0 - y1, r)
        : `M${x0},${y0}V${y1 - r}Q${x0},${y1} ${x0 + r},${y1}H${x0 + bw - r}Q${x0 + bw},${y1} ${x0 + bw},${y1 - r}V${y0}Z`; };
    const color = (d) => css(d.best.call === "CR" || d.best.call === "PR" ? "--good" : d.best.call === "SD" ? "--warning" : "--critical");
    svg.append("g").selectAll("path").data(data).join("path").attr("d", path).attr("fill", color);
    svg.append("line").attr("x1", m.l).attr("x2", W - m.r).attr("y1", y(0)).attr("y2", y(0)).attr("stroke", css("--axis"));
    svg.append("g").selectAll("rect").data(data).join("rect").attr("x", (d) => x(d.i)).attr("width", bw).attr("y", m.t).attr("height", H - m.t - m.b)
      .attr("fill", "transparent")
      .on("mousemove", (ev, d) => TIL.showTip(`<div class="tt-title">Virtual patient ${d.i}</div>
        <div class="row">Best response<b>${d.best.call} (${diamTxt(d.best.diam)})</b></div>
        <div class="row">Day 180<b>${d.final.call}${d.relapse ? ", relapsed" : ""}</b></div>
        <div class="row">Product stem-like share<b>${pct(d.stemFrac)}</b></div>
        <div class="row">Reactive clones<b>${pct(d.reactiveFrac)}</b></div>`, ev))
      .on("mouseleave", TIL.hideTip);
    drawScatter(data, color);
  }

  // Krishna inside the simulator: does product stemness predict the best response?
  function drawScatter(data, color) {
    const el = $("#scatter"); el.innerHTML = "";
    const W = el.clientWidth || 600, H = 240, m = { t: 10, r: 12, b: 36, l: 44 };
    const x = d3.scaleLinear().domain([0, d3.max(data, (d) => d.stemFrac) * 1.05]).nice().range([m.l, W - m.r]);
    const y = d3.scaleLinear().domain([-1, 1]).range([H - m.b, m.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Scatter of product stem-like share against best response");
    svg.append("g").attr("class", "gridline").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([-1, -0.3, 0.2, 1]).tickSize(-(W - m.l - m.r)).tickFormat(""));
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`)
      .call(d3.axisLeft(y).tickValues([-1, -0.3, 0, 0.2, 1]).tickFormat((d) => `${d > 0 ? "+" : ""}${Math.round(d * 100)}%`).tickSize(0).tickPadding(6)).call((g) => g.select(".domain").remove());
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`)
      .call(d3.axisBottom(x).ticks(6, "%").tickSize(0).tickPadding(7)).call((g) => g.select(".domain").remove());
    svg.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 4).attr("text-anchor", "middle").attr("font-size", 12)
      .attr("fill", css("--text-secondary")).text("Stem-like share of the patient's product");
    svg.append("g").selectAll("circle").data(data).join("circle").attr("cx", (d) => x(d.stemFrac)).attr("cy", (d) => y(d.v)).attr("r", 5.5)
      .attr("fill", color).attr("stroke", css("--surface-1")).attr("stroke-width", 2)
      .on("mousemove", (ev, d) => TIL.showTip(`<div class="tt-title">Virtual patient ${d.i}</div>
        <div class="row">Stem-like share<b>${pct(d.stemFrac)}</b></div><div class="row">Best response<b>${d.best.call} (${diamTxt(d.best.diam)})</b></div>
        <div class="row">Day 180<b>${d.final.call}${d.relapse ? ", relapsed" : ""}</b></div>`, ev))
      .on("mouseleave", TIL.hideTip);
    const rho = spearman(data.map((d) => d.stemFrac), data.map((d) => d.v));
    $("#scatterNote").textContent = `Spearman ρ = ${rho.toFixed(2)} between product stemness and best tumor change (negative = more stem-like, more shrinkage). ` +
      `Relapses: ${pct(d3.mean(data.filter((d) => d.stemFrac < d3.median(data, (e) => e.stemFrac)), (d) => +d.relapse))} of below-median-stemness products ` +
      `vs ${pct(d3.mean(data.filter((d) => d.stemFrac >= d3.median(data, (e) => e.stemFrac)), (d) => +d.relapse))} above.`;
  }
  function spearman(a, b) {
    const rank = (v) => { const o = v.map((x, i) => [x, i]).sort((p, q) => p[0] - q[0]), r = Array(v.length); o.forEach(([, i], k) => (r[i] = k)); return r; };
    const ra = rank(a), rb = rank(b), ma = d3.mean(ra), mb = d3.mean(rb);
    return d3.sum(ra, (x, i) => (x - ma) * (rb[i] - mb)) / Math.sqrt(d3.sum(ra, (x) => (x - ma) ** 2) * d3.sum(rb, (x) => (x - mb) ** 2));
  }

  let rt; window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(() => { if (last) { drawMuller(); drawStates(); drawTumor(); drawWaterfall(); } }, 150); });
  window.addEventListener("themechange", () => { if (last) { drawMuller(); drawStates(); drawTumor(); drawWaterfall(); } });
  run();
})();
