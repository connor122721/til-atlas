// Single-cell TILs: real datasets only. data/lab/index.json lists datasets exported by analysis/datasets/*/prepare.py.
(async function () {
  const { css, fmt, pct, esc } = TIL, { rankSumP, pTxt, strip } = TILStats;
  const $ = (s) => document.querySelector(s);
  const index = (await TIL.load("data/lab/index.json")).data;
  const catalog = (await TIL.load("data/catalog.json")).data;
  const receptors = (await TIL.load("data/lab/maynard2020/receptors.json")).data;
  const SERIES = ["--s1", "--s2", "--s3", "--s4", "--s5", "--s7"];
  const RESP = { Responder: "--s1", "Non-responder": "--s2" };
  const PHENO = { "stem-like": "--s1", effector: "--s3", terminal: "--s2" };
  const COMP_LABEL = { Tumor: "Tumor", NAT: "Normal adjacent", Blood: "Blood", tumor_pre: "Tumor (pre)", product: "Infusion product", blood_d30: "Blood, day 30", tumor_post: "Tumor (post)" };

  // ---------------- catalog: what else could be loaded
  const loaded = new Set((index?.datasets || []).map((d) => d.id).concat(receptors ? ["maynard2020"] : []));
  if (catalog) {
    const tag = { open: "good", check: "warning", unknown: "", none: "critical" };
    const q = catalog.datasets.filter((d) => d.use.includes("lab") && !loaded.has(d.id) && d.access !== "none");
    if (q.length) $("#queued").innerHTML = `<span class="muted">More public cohorts that can be added from your cluster:</span> ` + q.map((d) =>
      `<a href="datasets.html" class="pill ${tag[d.access] ?? ""}" title="${esc(d.status)}">${esc(d.short)}</a>`).join(" ") +
      ` <span class="muted">— the Datasets page says which already have a script.</span>`;
  }
  if (!index || !index.datasets.length) {
    ["#umap", "#strip"].forEach((s) => ($(s).innerHTML = `<p class="muted small">No single-cell dataset exported yet. Run <code>python analysis/datasets/wu2020/prepare.py</code>.</p>`));
    ["#expCard", "#cloneCard"].forEach((s) => ($(s).hidden = true));
  }

  let cur;
  const cache = {};
  async function loadDs(ds) {
    if (cache[ds.id]) return cache[ds.id];
    const get = (f) => TIL.load(`${ds.path}${f}.json`).then((r) => r.data);
    const [cells, pats, clones, exp] = await Promise.all([get("cells"), get("patients"), ds.hasClones ? get("clones") : null, ds.hasExpansion ? get("expansion") : null]);
    return (cache[ds.id] = { ds, cells, pats, clones, exp });
  }
  const groupOf = (d) => (cur.ds.groupBy === "cancer" ? d.group : d.response);
  function groupColors() {
    if (cur.ds.groupBy !== "cancer") return RESP;
    const gs = [...new Set(cur.pats.patients.map((p) => p.group))];
    return Object.fromEntries(gs.map((g, i) => [g, SERIES[i % SERIES.length]]));
  }

  async function select(id) {
    cur = await loadDs(index.datasets.find((d) => d.id === id));
    const { ds, cells } = cur;
    $("#dsNote").textContent = [ds.study, ds.therapy, `${cells.patients.length} patients`].filter(Boolean).join(" · ");
    $("#umapTitle").textContent = `CD8⁺ T cells · ${cells.meta.sample}`;
    $("#umapSub").textContent = `${fmt(cells.x.length)} cells drawn (up to 1,800 per patient) · UMAP and clusters from the original authors`;
    $("#clusterPick").innerHTML = cells.meta.clusters.map((c, i) => `<option value="${i}">${esc(c)}</option>`).join("");
    $("#colorBy option[value=compartment]").disabled = !cells.compartment;
    $("#metric option[value=dualExpanded]").disabled = !cur.exp;
    $("#cmpTitle").textContent = ds.groupBy === "cancer" ? "Patients by cancer type" : "Responders vs non-responders";
    $("#expCard").hidden = !cur.exp; $("#cloneCard").hidden = !cur.clones;
    if (cur.clones) $("#patient").innerHTML = Object.keys(cur.clones.patients).map((p) =>
      `<option value="${p}">${p} — ${cur.clones.patients[p].group || cur.clones.patients[p].response}</option>`).join("");
    all();
  }

  // ---------------- UMAP (canvas + quadtree hover)
  function drawUmap() {
    const cells = cur.cells, n = cells.x.length, el = $("#umap"); el.innerHTML = "";
    const W = el.clientWidth, H = el.clientHeight || W / 1.45, dpr = devicePixelRatio || 1, pad = 10;
    const cv = Object.assign(document.createElement("canvas"), { width: W * dpr, height: H * dpr });
    cv.style.width = W + "px"; cv.style.height = H + "px"; cv.setAttribute("role", "img"); cv.setAttribute("aria-label", "UMAP of CD8 T cells");
    el.append(cv);
    const ctx = cv.getContext("2d"); ctx.scale(dpr, dpr);
    const x = d3.scaleLinear().domain(d3.extent(cells.x)).range([pad, W - pad]);
    const y = d3.scaleLinear().domain(d3.extent(cells.y)).range([H - pad, pad]);
    const mode = $("#colorBy").value, pick = +$("#clusterPick").value, muted = css("--axis");
    let color, order = d3.range(n), legend = "";
    if (mode === "stem" || mode === "term") {
      const interp = mode === "stem" ? d3.interpolateRgb(css("--q100"), css("--q700")) : d3.interpolateRgb("#fbe0d5", "#7a2a0c");
      const v = cells[mode], lo = d3.quantile([...v].sort(d3.ascending), 0.02), hi = d3.quantile([...v].sort(d3.ascending), 0.98);
      const s = d3.scaleSequential(interp).domain([lo, hi]).clamp(true);
      color = (i) => s(v[i]); order.sort((a, b) => v[a] - v[b]);
      legend = `<span>low</span><span style="width:120px;height:8px;border-radius:4px;background:linear-gradient(90deg,${d3.range(0, 1.01, 0.25).map(interp).join(",")})"></span><span>high ${mode === "stem" ? "stem-like" : "terminal"} score</span>`;
    } else if (mode === "compartment") {
      const comps = cells.meta.compartmentsOfCells, cols = ["--s2", "--s1", "--s3"];
      color = (i) => css(cols[cells.compartment[i]]); d3.shuffle(order);
      legend = comps.map((c, j) => `<span><i style="background:${css(cols[j])}"></i>${esc(COMP_LABEL[c] || c)} · ${fmt(cells.compartment.filter((k) => k === j).length)}</span>`).join("");
    } else {
      color = (i) => (cells.cluster[i] === pick ? css("--s1") : muted);
      order.sort((a, b) => (cells.cluster[a] === pick) - (cells.cluster[b] === pick));
      const k = cells.cluster.filter((c) => c === pick).length;
      legend = `<span><i style="background:${css("--s1")}"></i>${esc(cells.meta.clusters[pick])} · ${fmt(k)} cells (${pct(k / n)})</span><span><i style="background:${muted}"></i>other CD8⁺ cells</span>`;
    }
    $("#umapLegend").innerHTML = legend;
    const r = W < 500 ? 1.4 : 1.9;
    order.forEach((i) => { ctx.fillStyle = color(i); ctx.beginPath(); ctx.arc(x(cells.x[i]), y(cells.y[i]), r, 0, 6.283); ctx.fill(); });
    const qt = d3.quadtree().x((i) => x(cells.x[i])).y((i) => y(cells.y[i])).addAll(d3.range(n));
    cv.addEventListener("mousemove", (ev) => {
      const b = cv.getBoundingClientRect(), i = qt.find(ev.clientX - b.left, ev.clientY - b.top, 10);
      if (i == null) return TIL.hideTip();
      const p = cells.patient[i];
      TIL.showTip(`<div class="tt-title">${esc(cells.meta.clusters[cells.cluster[i]])}</div>
        <div class="row">Patient<b>${p} · ${esc(cells.groups?.[p] || cells.responses?.[p] || "")}</b></div>
        ${cells.compartment ? `<div class="row">Tissue<b>${esc(COMP_LABEL[cells.meta.compartmentsOfCells[cells.compartment[i]]])}</b></div>` : ""}
        <div class="row">Stem-like score<b>${cells.stem[i].toFixed(2)}</b></div><div class="row">Terminal score<b>${cells.term[i].toFixed(2)}</b></div>`, ev);
    });
    cv.addEventListener("mouseleave", TIL.hideTip);
  }
  $("#colorBy").addEventListener("change", () => { $("#clusterPick").hidden = $("#colorBy").value !== "cluster"; drawUmap(); });
  $("#clusterPick").addEventListener("change", drawUmap);

  // ---------------- patient comparison
  function drawStrip() {
    const key = $("#metric").value, rows = cur.pats.patients, cols = groupColors();
    const groups = Object.keys(cols).map((g) => ({ key: g, label: g, color: cols[g], rows: rows.filter((d) => groupOf(d) === g) })).filter((g) => g.rows.length);
    strip($("#strip"), groups, (d) => d[key], { H: 280, fmtY: d3.format(".0%"),
      tip: (d) => `<div class="tt-title">${esc(d.patient)} · ${esc(groupOf(d))}</div><div class="row">Stem-like share<b>${pct(d.stemFrac, 1)}</b></div>
        <div class="row">Terminal share<b>${pct(d.termFrac, 1)}</b></div>${d.dualExpanded != null ? `<div class="row">In dual-expanded clones<b>${pct(d.dualExpanded, 1)}</b></div>` : ""}
        <div class="row">Tumor CD8⁺ cells<b>${fmt(d.nCells)}</b></div>` });
    if (groups.length === 2) {
      const a = groups[0].rows.map((d) => d[key]), b = groups[1].rows.map((d) => d[key]);
      $("#stats").innerHTML = `Two-sided Wilcoxon rank-sum, ${a.length} vs ${b.length} patients: <b>p = ${pTxt(rankSumP(a, b))}</b>.`;
    } else {
      $("#stats").textContent = "Medians per group: " + groups.map((g) => `${g.label} ${pct(d3.median(g.rows, (d) => d[key]))} (n=${g.rows.length})`).join(" · ") +
        ". Too few patients per cancer type for a test.";
    }
  }
  $("#metric").addEventListener("change", drawStrip);

  // ---------------- expansion categories (Wu)
  function drawExpansion() {
    if (!cur.exp) return;
    const el = $("#expansion"); el.innerHTML = "";
    const cats = [["Expanded in tumor and NAT", "--s1"], ["Expanded in tumor only", "--s2"], ["Tumor singleton", "--axis"]];
    const rows = cur.exp.patients;
    const W = el.clientWidth || 800, rowH = 22, lab = 120, right = 110, m = { t: 6, b: 24 };
    const H = m.t + rows.length * rowH + m.b;
    const x = d3.scaleLinear().domain([0, 1]).range([lab, W - right]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Share of tumor CD8 T cells by clone expansion pattern, per patient");
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(5, "%").tickSize(0).tickPadding(6)).call((g) => g.select(".domain").remove());
    rows.forEach((r, i) => {
      const y0 = m.t + i * rowH;
      svg.append("text").attr("x", lab - 8).attr("y", y0 + rowH / 2).attr("dy", "0.35em").attr("text-anchor", "end").attr("font-size", 12)
        .attr("fill", css("--text-secondary")).text(`${r.patient} · ${r.group}`);
      let acc = 0;
      cats.forEach(([c, col]) => {
        const v = r.categories[c] || 0;
        svg.append("rect").attr("x", x(acc)).attr("y", y0 + 3).attr("width", Math.max(0, x(acc + v) - x(acc) - 2)).attr("height", rowH - 6)
          .attr("fill", css(col)).attr("rx", 2)
          .on("mousemove", (ev) => TIL.showTip(`<div class="tt-title">${esc(r.patient)} · ${esc(r.group)}</div>` + cats.map(([cc, co]) =>
            `<div class="row"><i style="background:${css(co)}"></i>${cc}<b>${pct(r.categories[cc], 1)}</b></div>`).join("") +
            `<div class="row">Tumor CD8⁺ cells with TCR<b>${fmt(r.tumorCells)}</b></div>` +
            (r.hasBlood ? `<div class="row">In clones also seen in blood<b>${pct(r.sharedWithBlood, 1)}</b></div>` : ""), ev))
          .on("mouseleave", TIL.hideTip);
        acc += v;
      });
      if (r.hasBlood) svg.append("text").attr("x", W - right + 8).attr("y", y0 + rowH / 2).attr("dy", "0.35em").attr("font-size", 11.5)
        .attr("fill", css("--text-secondary")).text(`${pct(r.sharedWithBlood)} also in blood`);
    });
    $("#expLegend").innerHTML = cats.map(([c, col]) => `<span><i style="background:${css(col)}"></i>${c}</span>`).join("") +
      `<span class="muted">· expanded = 2 or more cells; NAT = normal adjacent tissue</span>`;
  }

  // ---------------- clone slope chart
  function drawSlope() {
    if (!cur.clones) return;
    const clones = cur.clones, el = $("#slope"); el.innerHTML = "";
    const pid = $("#patient").value, P = clones.patients[pid];
    const comps = clones.meta.compartments.filter((c) => P.clones.some((k) => c in k.freq));
    $("#patientNote").textContent = `${P.clones.length} clones · ${P.clones.filter((c) => c.reactive).length} predicted tumor-reactive`;
    $("#cloneLegend").innerHTML = Object.entries(PHENO).map(([k, c]) => `<span><i style="background:${css(c)}"></i>Mostly ${k} in tumor</span>`).join("") +
      `<span><i style="background:transparent;border:2px solid ${css("--text-secondary")}"></i>predicted tumor-reactive (thick line)</span>`;
    const W = el.clientWidth || 800, H = 340, m = { t: 14, r: 24, b: 30, l: 54 };
    const x = d3.scalePoint().domain(comps).range([m.l, W - m.r]).padding(0.15);
    const floor = -4.3, y = d3.scaleLinear().domain([floor, 0]).range([H - m.b, m.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Clone frequencies across tissues");
    svg.append("rect").attr("x", m.l).attr("width", W - m.l - m.r).attr("y", y(-4)).attr("height", y(floor) - y(-4)).attr("fill", css("--surface-2"));
    svg.append("text").attr("x", m.l + 4).attr("y", (y(-4) + y(floor)) / 2).attr("dy", "0.35em").attr("font-size", 10.5).attr("fill", css("--text-muted")).text("not detected");
    const ticks = [-4, -3, -2, -1, 0];
    svg.append("g").attr("class", "gridline").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues(ticks).tickSize(-(W - m.l - m.r)).tickFormat(""));
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`)
      .call(d3.axisLeft(y).tickValues(ticks).tickFormat((d) => ({ 0: "100%", "-1": "10%", "-2": "1%", "-3": "0.1%", "-4": "0.01%" }[d])).tickSize(0).tickPadding(6))
      .call((g) => g.select(".domain").remove());
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickFormat((d) => COMP_LABEL[d] || d).tickSize(0).tickPadding(10))
      .call((g) => g.select(".domain").remove());
    const yv = (v) => y(v == null ? (floor - 4) / 2 : Math.max(v, -4));
    const line = d3.line().x((d) => x(d[0])).y((d) => yv(d[1]));
    const paths = svg.append("g").selectAll("g").data(P.clones).join("g");
    paths.append("path").attr("d", (c) => line(comps.map((k) => [k, c.freq[k]]))).attr("fill", "none")
      .attr("stroke", (c) => css(PHENO[c.phenotype])).attr("stroke-width", (c) => (c.reactive ? 2.6 : 1.2)).attr("stroke-opacity", 0.85);
    paths.selectAll("circle").data((c) => comps.map((k) => ({ c, k, v: c.freq[k] }))).join("circle")
      .attr("cx", (d) => x(d.k)).attr("cy", (d) => yv(d.v)).attr("r", 3.5).attr("fill", (d) => css(PHENO[d.c.phenotype]))
      .attr("stroke", css("--surface-1")).attr("stroke-width", 1.5);
    paths.append("path").attr("d", (c) => line(comps.map((k) => [k, c.freq[k]]))).attr("fill", "none").attr("stroke", "transparent").attr("stroke-width", 12)
      .on("mousemove", (ev, c) => {
        paths.style("opacity", (o) => (o === c ? 1 : 0.15));
        TIL.showTip(`<div class="tt-title mono">${esc(c.id || c.cdr3)}</div><div class="row">${c.reactive ? "Predicted tumor-reactive" : "Not flagged"} · mostly ${c.phenotype}</div>` +
          comps.map((k) => `<div class="row">${COMP_LABEL[k] || k}<b>${c.freq[k] == null ? "not detected" : pct(10 ** c.freq[k], 2)}</b></div>`).join(""), ev);
      })
      .on("mouseleave", () => { paths.style("opacity", 1); TIL.hideTip(); });
  }
  $("#patient").addEventListener("change", drawSlope);

  // ---------------- BCR / TCR panel (Maynard 2020)
  function drawReceptors() {
    if (!receptors) { $("#bcrCard").hidden = true; return; }
    const R = receptors, ST = R.meta.stages, stages = Object.keys(ST);
    $("#bcrSub").textContent = `${R.meta.study} · ${fmt(R.meta.nBCR)} B cell and ${fmt(R.meta.nTCR)} T cell receptors from lung cancer biopsies taken before and during ${R.meta.therapy}.`;
    const ISO = [["IgG", "--s1"], ["IgA", "--s3"], ["IgM/IgD", "--s4"]];
    const el = $("#iso"); el.innerHTML = "";
    const W = el.clientWidth || 380, rowH = 34, lab = 150, H = stages.length * rowH + 26;
    const x = d3.scaleLinear().domain([0, 1]).range([lab, W - 10]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Antibody class by treatment stage");
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - 22})`).call(d3.axisBottom(x).ticks(4, "%").tickSize(0).tickPadding(6)).call((g) => g.select(".domain").remove());
    stages.forEach((st, i) => {
      const p = R.pooled[st], tot = d3.sum(Object.values(p.isotype)), y0 = i * rowH;
      svg.append("text").attr("x", lab - 8).attr("y", y0 + rowH / 2 - 5).attr("text-anchor", "end").attr("font-size", 12).attr("fill", css("--text-primary")).text(ST[st].replace(/ \(.*/, ""));
      svg.append("text").attr("x", lab - 8).attr("y", y0 + rowH / 2 + 9).attr("text-anchor", "end").attr("font-size", 11).attr("fill", css("--text-muted")).text(`${fmt(tot)} B cells`);
      let acc = 0;
      ISO.forEach(([k, col]) => {
        const v = (p.isotype[k] || 0) / tot;
        svg.append("rect").attr("x", x(acc)).attr("y", y0 + 6).attr("width", Math.max(0, x(acc + v) - x(acc) - 2)).attr("height", rowH - 12).attr("rx", 2).attr("fill", css(col))
          .on("mousemove", (ev) => TIL.showTip(`<div class="tt-title">${esc(ST[st])}</div>` + ISO.map(([kk, cc]) =>
            `<div class="row"><i style="background:${css(cc)}"></i>${kk}<b>${pct((p.isotype[kk] || 0) / tot, 1)} (${p.isotype[kk] || 0})</b></div>`).join(""), ev))
          .on("mouseleave", TIL.hideTip);
        if (v > 0.12) svg.append("text").attr("x", x(acc) + 6).attr("y", y0 + rowH / 2).attr("dy", "0.35em").attr("font-size", 11).attr("fill", col === "--s4" ? "#0b0b0b" : "#fff").attr("pointer-events", "none").text(pct(v));
        acc += v;
      });
    });
    $("#isoLegend").innerHTML = ISO.map(([k, c]) => `<span><i style="background:${css(c)}"></i>${k}</span>`).join("");

    const metric = $("#bcrMetric").value;
    const val = { igg: (s) => (s.bcr.n >= 10 ? s.bcr.isotype?.IgG : null), bexp: (s) => (s.bcr.n >= 10 ? s.bcr.expandedShare : null), texp: (s) => (s.tcr.n >= 10 ? s.tcr.expandedShare : null) }[metric];
    const groups = stages.map((st, i) => ({ key: st, label: ST[st].replace(/ \(.*/, "").replace("Residual disease on therapy", "Residual disease"), color: ["--s1", "--s3", "--s2"][i],
      rows: R.samples.filter((s) => s.stage === st && val(s) != null) }));
    strip($("#bcrStrip"), groups, val, { H: 230, fmtY: d3.format(".0%"),
      tip: (s) => `<div class="tt-title">${esc(s.sample)} · ${esc(s.patient)}</div><div class="row">${esc(s.tissue)}${s.response ? ` · best response ${esc(s.response)}` : ""}</div>
        <div class="row">B cells with BCR<b>${s.bcr.n}</b></div><div class="row">T cells with TCR<b>${s.tcr.n}</b></div>
        ${s.bcr.isotype ? `<div class="row">IgG / IgA / IgM-D<b>${pct(s.bcr.isotype.IgG)} / ${pct(s.bcr.isotype.IgA)} / ${pct(s.bcr.isotype["IgM/IgD"])}</b></div>` : ""}` });
    $("#bcrStats").textContent = "One dot per biopsy with ≥10 receptors: " + groups.map((g) => `${g.label} ${g.rows.length}`).join(", ") +
      ". Several biopsies can come from one patient, so this is descriptive.";
  }
  $("#bcrMetric").addEventListener("change", drawReceptors);

  const all = () => { if (cur) { drawUmap(); drawStrip(); drawExpansion(); drawSlope(); } drawReceptors(); };
  let t; window.addEventListener("resize", () => { clearTimeout(t); t = setTimeout(all, 150); });
  window.addEventListener("themechange", all);
  if (index?.datasets.length) {
    index.datasets.forEach((d) => $("#dataset").insertAdjacentHTML("beforeend", `<option value="${esc(d.id)}">${esc(d.label)}</option>`));
    select(index.datasets[0].id);
  } else drawReceptors();
})();
