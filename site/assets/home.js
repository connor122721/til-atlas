// Home: real numbers and a real result (CD8 signal vs checkpoint response across cohorts).
(async function () {
  const { css, fmt, esc } = TIL;
  const $ = (s) => document.querySelector(s);
  const [trials, bulkIdx, bulkSum, lab] = await Promise.all(["data/trials.json", "data/bulk/index.json", "data/bulk/summary.json", "data/lab/wu2020/patients.json"]
    .map((p) => TIL.load(p).then((r) => r.data)));

  const tiles = [];
  if (trials) tiles.push(["Trials tracked", fmt(trials.trials.length), `ClinicalTrials.gov, ${trials.meta.generated}`]);
  if (bulkIdx) tiles.push(["Patients with tumor RNA-seq", fmt(d3.sum(bulkIdx.cohorts, (c) => c.n)), `${bulkIdx.cohorts.length} immunotherapy cohorts`]);
  if (lab) tiles.push(["Tumor CD8⁺ T cells profiled", fmt(d3.sum(lab.patients, (p) => p.nCells)), `${lab.patients.length} patients, Wu et al. 2020`]);
  $("#heroStats").style.gridTemplateColumns = `repeat(${Math.min(3, tiles.length)}, minmax(0, 1fr))`;
  $("#heroStats").innerHTML = tiles.slice(0, 3).map(([l, v, n]) => `<div class="card stat"><div class="label">${l}</div><div class="value">${v}</div><div class="note">${n}</div></div>`).join("");

  function art() {
    const el = $("#heroArt"); el.innerHTML = "";
    if (!bulkIdx || !bulkSum) return;
    const rows = bulkIdx.cohorts.map((c) => ({ c, s: bulkSum.stats[c.id]?.response?.cd8 })).filter((r) => r.s);
    const W = el.clientWidth || 460, rowH = 30, lab = 112, m = { t: 26, r: 14, b: 28 }, H = m.t + rows.length * rowH + m.b;
    const x = d3.scaleLinear().domain([0.2, 1]).range([lab, W - m.r]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img")
      .attr("aria-label", "AUC of CD8 T cell expression for objective response, by immunotherapy cohort");
    svg.append("g").attr("class", "gridline").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([0.25, 0.5, 0.75, 1]).tickSize(-(H - m.t - m.b)).tickFormat(""));
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([0.25, 0.5, 0.75, 1]).tickFormat(d3.format(".2~f")).tickSize(0).tickPadding(7)).call((g) => g.select(".domain").remove());
    svg.append("line").attr("x1", x(0.5)).attr("x2", x(0.5)).attr("y1", m.t - 6).attr("y2", H - m.b).attr("stroke", css("--text-muted")).attr("stroke-dasharray", "4 3");
    svg.append("text").attr("x", x(0.5) + 6).attr("y", m.t - 10).attr("font-size", 11).attr("fill", css("--text-muted")).text("higher in responders →");
    const r = d3.scaleSqrt().domain([0, d3.max(rows, (d) => d.s.nR + d.s.nNR)]).range([3, 8]);
    rows.forEach((d, i) => {
      const y = m.t + i * rowH + rowH / 2;
      svg.append("text").attr("x", 0).attr("y", y).attr("dy", "-0.1em").attr("font-size", 12).attr("font-weight", 600).attr("fill", css("--text-primary")).text(d.c.label);
      svg.append("text").attr("x", 0).attr("y", y).attr("dy", "1.05em").attr("font-size", 10.5).attr("fill", css("--text-muted")).text(d.c.cancer);
      svg.append("line").attr("x1", x(d.s.lo)).attr("x2", x(d.s.hi)).attr("y1", y).attr("y2", y).attr("stroke", css("--s1")).attr("stroke-width", 2);
      svg.append("circle").attr("cx", x(d.s.auc)).attr("cy", y).attr("r", r(d.s.nR + d.s.nNR)).attr("fill", css("--s1")).attr("stroke", css("--surface-1")).attr("stroke-width", 2)
        .on("mousemove", (ev) => TIL.showTip(`<div class="tt-title">${esc(d.c.label)} · ${esc(d.c.cancer)}</div><div class="row">AUC<b>${d.s.auc.toFixed(2)} (${d.s.lo.toFixed(2)}–${d.s.hi.toFixed(2)})</b></div><div class="row">Patients<b>${d.s.nR} R / ${d.s.nNR} NR</b></div>`, ev))
        .on("mouseleave", TIL.hideTip);
    });
  }

  function flow() {
    const el = $("#flow"); el.innerHTML = "";
    const steps = [["Weekly cron", "GitHub Action, Mondays"], ["Query APIs", "ClinicalTrials.gov + GEO"], ["Normalize + tag", "regex rules in tags.json"],
      ["Diff vs last week", "new · status changes · removed"], ["Commit JSON", "site/data/*.json"], ["Pages deploy", "static site rebuilds"]];
    const W = el.clientWidth, narrow = W < 640, cols = narrow ? 2 : 6, bw = (W - (cols - 1) * 22) / cols, bh = 62, rowsN = Math.ceil(steps.length / cols);
    const H = rowsN * bh + (rowsN - 1) * 24 + 4;
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Pipeline: cron, API, normalize, diff, commit, deploy");
    steps.forEach(([a, b], i) => {
      const cx = (i % cols) * (bw + 22), cy = Math.floor(i / cols) * (bh + 24) + 2;
      svg.append("rect").attr("x", cx).attr("y", cy).attr("width", bw).attr("height", bh).attr("rx", 8).attr("fill", css("--surface-2")).attr("stroke", css("--border"));
      svg.append("text").attr("x", cx + 12).attr("y", cy + 25).attr("font-size", 13.5).attr("font-weight", 650).attr("fill", css("--text-primary")).text(a);
      svg.append("text").attr("x", cx + 12).attr("y", cy + 45).attr("font-size", 11.5).attr("fill", css("--text-secondary")).text(b);
      if (i % cols < cols - 1 && i < steps.length - 1)
        svg.append("path").attr("d", `M${cx + bw + 4},${cy + bh / 2}h14m-5,-4l5,4l-5,4`).attr("fill", "none").attr("stroke", css("--text-muted")).attr("stroke-width", 1.5);
    });
  }
  const all = () => { art(); flow(); };
  let t; window.addEventListener("resize", () => { clearTimeout(t); t = setTimeout(all, 150); });
  window.addEventListener("themechange", all);
  all();
})();
