// Replication lab: pick a dataset (data/lab/index.json), explore its CD8 cells, compare responders,
// compare every loaded dataset side by side, and follow clones where TCR data exist.
(async function () {
  const { css, fmt, pct, esc } = TIL;
  const $ = (s) => document.querySelector(s);
  const index = (await TIL.load("data/lab/index.json")).data;
  const catalog = (await TIL.load("data/catalog.json")).data;
  if (!index || !index.datasets.length) {
    TIL.banner($("#main"), "No lab data found.", "Run <code>python analysis/make_synthetic_demo.py site/data/lab</code> or the R exports."); return;
  }
  const RESP = { Responder: "--s1", "Non-responder": "--s2" };
  const cache = {};
  async function loadDs(ds) {
    if (cache[ds.id]) return cache[ds.id];
    const [c, p, k] = await Promise.all(["cells", "patients", ds.hasClones ? "clones" : null]
      .map((f) => (f ? TIL.load(`${ds.path}${f}.json`).then((r) => r.data) : null)));
    return (cache[ds.id] = { ds, cells: c, pats: p, clones: k });
  }
  // patient tables for every dataset (small) feed the cross-dataset panel
  const allPats = await Promise.all(index.datasets.map((ds) => TIL.load(`${ds.path}patients.json`).then((r) => ({ ds, pats: r.data }))));

  // ---------------- dataset picker + what's queued from the catalog
  index.datasets.forEach((d) => $("#dataset").insertAdjacentHTML("beforeend",
    `<option value="${esc(d.id)}">${esc(d.label)}${d.synthetic ? " (synthetic)" : ""}</option>`));
  const loadedIds = new Set(index.datasets.map((d) => d.id));
  if (catalog) {
    const queued = catalog.datasets.filter((d) => d.use.includes("lab") && !loadedIds.has(d.id));
    const tag = { open: "good", check: "warning", unknown: "" };
    $("#queued").innerHTML = queued.length ? `<span class="muted">Not loaded yet:</span> ` + queued.map((d) =>
      `<a href="datasets.html" class="pill ${tag[d.access]}" title="${esc(d.status)}">${esc(d.short)}</a>`).join(" ") +
      ` <span class="muted">— run their pipelines in <code>analysis/datasets/</code> to add them here.</span>` : "";
  }

  let cur;
  async function select(id) {
    cur = await loadDs(index.datasets.find((d) => d.id === id));
    const { ds, cells } = cur;
    $("#dsBanner").innerHTML = "";
    if (ds.synthetic) TIL.banner($("#dsBanner"), "Synthetic demo data.",
      "Every cell, patient and clone in this dataset is invented (the responder difference was put in by hand) so you can see how the views work. " +
      "Pick a real dataset above once its pipeline has run.");
    $("#dsNote").textContent = [ds.study, ds.therapy !== "—" ? ds.therapy : null, `${cells.patients.length} patients`, ds.sample].filter(Boolean).join(" · ");
    $("#umapTitle").textContent = `CD8⁺ T cells · ${ds.sample}`;
    $("#umapSub").textContent = `${fmt(cells.x.length)} CD8⁺ cells from ${cells.patients.length} samples · ${cells.meta.sample}`;
    $("#clusterPick").innerHTML = cells.meta.clusters.map((c, i) => `<option value="${i}">${esc(c)}</option>`).join("");
    const hasPers = cur.pats.patients.some((p) => p.persistence30 != null);
    $("#metric option[value=persistence30]").disabled = !hasPers;
    if (!hasPers && $("#metric").value === "persistence30") $("#metric").value = "stemFrac";
    $("#cloneCard").hidden = !cur.clones;
    if (cur.clones) {
      $("#patient").innerHTML = Object.keys(cur.clones.patients).map((p) =>
        `<option value="${p}">${p} — ${cur.clones.patients[p].response}</option>`).join("");
    }
    all();
  }
  $("#dataset").addEventListener("change", (e) => select(e.target.value));

  // ---------------- UMAP on canvas (thousands of points), hover via quadtree
  function drawUmap() {
    const cells = cur.cells, n = cells.x.length, el = $("#umap"); el.innerHTML = "";
    const W = el.clientWidth, H = el.clientHeight || W / 1.45, dpr = devicePixelRatio || 1, pad = 10;
    const cv = Object.assign(document.createElement("canvas"), { width: W * dpr, height: H * dpr });
    cv.style.width = W + "px"; cv.style.height = H + "px"; cv.setAttribute("role", "img"); cv.setAttribute("aria-label", "UMAP scatter of CD8 T cells");
    el.append(cv);
    const ctx = cv.getContext("2d"); ctx.scale(dpr, dpr);
    const x = d3.scaleLinear().domain(d3.extent(cells.x)).range([pad, W - pad]);
    const y = d3.scaleLinear().domain(d3.extent(cells.y)).range([H - pad, pad]);
    const mode = $("#colorBy").value, pick = +$("#clusterPick").value, muted = css("--axis");
    let color, order = d3.range(n), legend = "";
    if (mode === "stem" || mode === "term") {
      const interp = mode === "stem" ? d3.interpolateRgb(css("--q100"), css("--q700")) : d3.interpolateRgb("#fbe0d5", "#7a2a0c");
      const s = d3.scaleSequential(interp).domain(d3.extent(cells[mode]));
      color = (i) => s(cells[mode][i]);
      order.sort((a, b) => cells[mode][a] - cells[mode][b]);
      const stops = d3.range(0, 1.01, 0.25).map((t) => interp(t)).join(",");
      legend = `<span>low</span><span style="width:120px;height:8px;border-radius:4px;background:linear-gradient(90deg,${stops})"></span><span>high ${mode === "stem" ? "stem-like" : "terminal"} score (UCell)</span>`;
    } else if (mode === "response") {
      color = (i) => css(RESP[cells.responses[cells.patient[i]]]);
      d3.shuffle(order);
      legend = Object.entries(RESP).map(([k, c]) => `<span><i style="background:${css(c)}"></i>${k}</span>`).join("");
    } else {
      color = (i) => (cells.cluster[i] === pick ? css("--s1") : muted);
      order.sort((a, b) => (cells.cluster[a] === pick) - (cells.cluster[b] === pick));
      const k = cells.cluster.filter((c) => c === pick).length;
      legend = `<span><i style="background:${css("--s1")}"></i>${esc(cells.meta.clusters[pick])} · ${fmt(k)} cells (${pct(k / n)})</span><span><i style="background:${muted}"></i>other cells</span>`;
    }
    $("#umapLegend").innerHTML = legend;
    const r = W < 500 ? 1.6 : 2.2;
    order.forEach((i) => { ctx.fillStyle = color(i); ctx.beginPath(); ctx.arc(x(cells.x[i]), y(cells.y[i]), r, 0, 6.283); ctx.fill(); });
    const qt = d3.quadtree().x((i) => x(cells.x[i])).y((i) => y(cells.y[i])).addAll(d3.range(n));
    cv.addEventListener("mousemove", (ev) => {
      const b = cv.getBoundingClientRect(), i = qt.find(ev.clientX - b.left, ev.clientY - b.top, 10);
      if (i == null) return TIL.hideTip();
      const p = cells.patient[i];
      TIL.showTip(`<div class="tt-title">${esc(cells.meta.clusters[cells.cluster[i]])}</div>
        <div class="row">Patient<b>${p} · ${cells.responses[p]}</b></div>
        <div class="row">Stem-like score<b>${cells.stem[i].toFixed(2)}</b></div>
        <div class="row">Terminal score<b>${cells.term[i].toFixed(2)}</b></div>`, ev);
    });
    cv.addEventListener("mouseleave", TIL.hideTip);
  }
  $("#colorBy").addEventListener("change", () => { $("#clusterPick").hidden = $("#colorBy").value !== "cluster"; drawUmap(); });
  $("#clusterPick").addEventListener("change", drawUmap);

  // ---------------- Wilcoxon rank-sum: exact by enumeration for small n, normal approximation otherwise
  function rankSumP(a, b) {
    const all = [...a, ...b], ranks = (() => { const s = all.map((v, i) => [v, i]).sort((p, q) => p[0] - q[0]), r = Array(all.length);
      for (let i = 0; i < s.length;) { let j = i; while (j + 1 < s.length && s[j + 1][0] === s[i][0]) j++;
        for (let k = i; k <= j; k++) r[s[k][1]] = (i + j) / 2 + 1; i = j + 1; } return r; })();
    const obs = d3.sum(ranks.slice(0, a.length)), N = all.length, mean = a.length * (N + 1) / 2;
    if (!a.length || !b.length) return NaN;
    if (N > 22) {
      const sd = Math.sqrt(a.length * b.length * (N + 1) / 12), z = (Math.abs(obs - mean) - 0.5) / sd;
      return Math.min(1, 2 * (1 - normCdf(z)));
    }
    let extreme = 0, total = 0;
    (function rec(start, k, sum) {
      if (k === 0) { total++; if (Math.abs(sum - mean) >= Math.abs(obs - mean) - 1e-9) extreme++; return; }
      for (let i = start; i <= N - k; i++) rec(i + 1, k - 1, sum + ranks[i]);
    })(0, a.length, 0);
    return extreme / total;
  }
  const normCdf = (z) => 0.5 * (1 + Math.tanh(Math.sqrt(2 / Math.PI) * (z + 0.044715 * z ** 3)));
  const pTxt = (p) => (isNaN(p) ? "n/a" : p < 0.001 ? p.toExponential(1) : p.toFixed(3));
  const tip = (d, ds) => `<div class="tt-title">${esc(d.patient)} · ${d.response}</div>${ds ? `<div class="row">${esc(ds.label)}</div>` : ""}
    <div class="row">Stem-like share<b>${pct(d.stemFrac, 1)}</b></div><div class="row">Terminal share<b>${pct(d.termFrac, 1)}</b></div>
    ${d.persistence30 != null ? `<div class="row">Persistence d30<b>${pct(d.persistence30)}</b></div>` : ""}<div class="row">CD8⁺ cells<b>${fmt(d.nCells)}</b></div>`;

  /** Strip plot of patients by response. groups: [{key, label, rows}] with rows = patient records. */
  function strip(el, groups, key, { H = 280, labelMedians = true, dsOf = null } = {}) {
    el.innerHTML = "";
    const rowsAll = groups.flatMap((g) => g.rows).filter((d) => d[key] != null);
    if (!rowsAll.length) { el.innerHTML = `<p class="muted small">This metric isn't available for these datasets.</p>`; return; }
    const W = el.clientWidth || 300, m = { t: 10, r: 10, b: dsOf ? 44 : 28, l: 42 };
    const x = d3.scaleBand().domain(groups.map((g) => g.key)).range([m.l, W - m.r]).padding(0.3);
    const y = d3.scaleLinear().domain([0, Math.max(0.1, d3.max(rowsAll, (d) => d[key]) * 1.12)]).nice().range([H - m.b, m.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Per-patient values by response group");
    svg.append("g").attr("class", "gridline").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickSize(-(W - m.l - m.r)).tickFormat(""));
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5, "%").tickSize(0).tickPadding(6)).call((g) => g.select(".domain").remove());
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`)
      .call(d3.axisBottom(x).tickFormat((k) => groups.find((g) => g.key === k).label).tickSize(0).tickPadding(8)).call((g) => g.select(".domain").remove());
    if (dsOf) { // dataset names under each responder/non-responder pair
      d3.groups(groups, (g) => g.ds.id).forEach(([, gs]) => {
        const x0 = x(gs[0].key), x1 = x(gs[gs.length - 1].key) + x.bandwidth();
        svg.append("text").attr("x", (x0 + x1) / 2).attr("y", H - 6).attr("text-anchor", "middle").attr("font-size", 11.5)
          .attr("font-weight", 600).attr("fill", css("--text-secondary")).text(gs[0].ds.label);
      });
    }
    const jr = d3.randomUniform.source(d3.randomLcg(0.37))(-1, 1);
    groups.forEach((g) => {
      const vals = g.rows.map((d) => d[key]).filter((v) => v != null); if (!vals.length) return;
      const med = d3.median(vals);
      svg.append("line").attr("x1", x(g.key)).attr("x2", x(g.key) + x.bandwidth()).attr("y1", y(med)).attr("y2", y(med))
        .attr("stroke", css("--text-primary")).attr("stroke-width", 2);
      if (labelMedians) svg.append("text").attr("x", x(g.key) + x.bandwidth() + 4).attr("y", y(med)).attr("dy", "0.35em").attr("font-size", 11)
        .attr("fill", css("--text-secondary")).text(`median ${pct(med)}`);
      svg.append("g").selectAll("circle").data(g.rows.filter((d) => d[key] != null)).join("circle")
        .attr("cx", () => x(g.key) + x.bandwidth() / 2 + jr() * x.bandwidth() * 0.32).attr("cy", (d) => y(d[key])).attr("r", 5)
        .attr("fill", (d) => css(RESP[d.response])).attr("stroke", css("--surface-1")).attr("stroke-width", 2)
        .on("mousemove", (ev, d) => TIL.showTip(tip(d, g.ds), ev)).on("mouseleave", TIL.hideTip);
    });
  }
  const byResp = (rows, ds) => Object.keys(RESP).map((r) => ({ key: `${ds ? ds.id + "|" : ""}${r}`, label: r === "Responder" ? (ds ? "R" : "Responder") : (ds ? "NR" : "Non-responder"),
    rows: rows.filter((d) => d.response === r), ds }));

  function drawStrip() {
    const key = $("#metric").value, rows = cur.pats.patients;
    strip($("#strip"), byResp(rows), key);
    const a = rows.filter((d) => d.response === "Responder" && d[key] != null).map((d) => d[key]);
    const b = rows.filter((d) => d.response === "Non-responder" && d[key] != null).map((d) => d[key]);
    $("#stats").innerHTML = `Two-sided Wilcoxon rank-sum, ${a.length} vs ${b.length} patients: <b>p = ${pTxt(rankSumP(a, b))}</b>` +
      `${a.length + b.length <= 22 ? " (exact, by enumerating every relabeling)" : " (normal approximation)"}.`;
  }
  function drawPooled() {
    const key = $("#metric").value;
    const groups = allPats.filter((x) => x.pats).flatMap(({ ds, pats }) => byResp(pats.patients, ds));
    strip($("#pooled"), groups, key, { H: 300, labelMedians: false, dsOf: true });
    $("#pooledStats").innerHTML = allPats.filter((x) => x.pats).map(({ ds, pats }) => {
      const r = pats.patients.filter((d) => d[key] != null);
      const a = r.filter((d) => d.response === "Responder").map((d) => d[key]), b = r.filter((d) => d.response === "Non-responder").map((d) => d[key]);
      return `<b>${esc(ds.label)}</b>${ds.synthetic ? " (synthetic)" : ""}: ${a.length} vs ${b.length}, p = ${pTxt(rankSumP(a, b))}`;
    }).join(" · ");
  }
  $("#metric").addEventListener("change", () => { drawStrip(); drawPooled(); });

  // ---------------- clone slope chart
  const PHENO = { "stem-like": "--s1", terminal: "--s2" };
  const LABEL = { tumor_pre: "Tumor (pre)", product: "Infusion product", blood_d30: "Blood, day 30", tumor_post: "Tumor (post)" };
  function drawSlope() {
    if (!cur.clones) return;
    const clones = cur.clones, el = $("#slope"); el.innerHTML = "";
    const pid = $("#patient").value, P = clones.patients[pid], comps = clones.meta.compartments;
    $("#patientNote").textContent = `${P.clones.length} clonotypes · ${P.clones.filter((c) => c.reactive).length} flagged tumor-reactive`;
    $("#cloneLegend").innerHTML = Object.entries(PHENO).map(([k, c]) => `<span><i style="background:${css(c)}"></i>${k === "stem-like" ? "Mostly stem-like in product" : "Mostly terminal in product"}</span>`).join("") +
      `<span><i style="background:transparent;border:2px solid ${css("--text-secondary")}"></i>tumor-reactive (thicker line)</span>`;
    const W = el.clientWidth || 800, H = 340, m = { t: 14, r: 20, b: 30, l: 48 };
    const x = d3.scalePoint().domain(comps).range([m.l, W - m.r]).padding(0.1);
    const floor = -5.3, y = d3.scaleLinear().domain([floor, 0]).range([H - m.b, m.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Clone frequencies across samples");
    svg.append("rect").attr("x", m.l).attr("width", W - m.l - m.r).attr("y", y(-5)).attr("height", y(floor) - y(-5)).attr("fill", css("--surface-2"));
    svg.append("text").attr("x", m.l + 4).attr("y", (y(-5) + y(floor)) / 2).attr("dy", "0.35em").attr("font-size", 10.5).attr("fill", css("--text-muted")).text("not detected");
    svg.append("g").attr("class", "gridline").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([-5, -4, -3, -2, -1, 0]).tickSize(-(W - m.l - m.r)).tickFormat(""));
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`)
      .call(d3.axisLeft(y).tickValues([-5, -4, -3, -2, -1, 0]).tickFormat((d) => ({ 0: "100%", "-1": "10%", "-2": "1%", "-3": "0.1%", "-4": "0.01%", "-5": "0.001%" }[d])).tickSize(0).tickPadding(6))
      .call((g) => g.select(".domain").remove());
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickFormat((d) => LABEL[d] || d).tickSize(0).tickPadding(10))
      .call((g) => g.select(".domain").remove());
    const yv = (v) => y(v == null ? (floor - 5) / 2 : v);
    const line = d3.line().x((d) => x(d[0])).y((d) => yv(d[1]));
    const paths = svg.append("g").selectAll("g").data(P.clones).join("g");
    paths.append("path").attr("d", (c) => line(comps.map((k) => [k, c.freq[k]]))).attr("fill", "none")
      .attr("stroke", (c) => css(PHENO[c.phenotype])).attr("stroke-width", (c) => (c.reactive ? 2.5 : 1.2)).attr("stroke-opacity", 0.85);
    paths.selectAll("circle").data((c) => comps.map((k) => ({ c, k, v: c.freq[k] }))).join("circle")
      .attr("cx", (d) => x(d.k)).attr("cy", (d) => yv(d.v)).attr("r", 3.5).attr("fill", (d) => css(PHENO[d.c.phenotype]))
      .attr("stroke", css("--surface-1")).attr("stroke-width", 1.5);
    paths.append("path").attr("d", (c) => line(comps.map((k) => [k, c.freq[k]]))).attr("fill", "none").attr("stroke", "transparent").attr("stroke-width", 12)
      .on("mousemove", function (ev, c) {
        paths.style("opacity", (o) => (o === c ? 1 : 0.15));
        TIL.showTip(`<div class="tt-title mono">${esc(c.cdr3)}</div>
          <div class="row">${c.reactive ? "Tumor-reactive" : "Reactivity not shown"} · ${c.phenotype}</div>` +
          comps.map((k) => `<div class="row">${LABEL[k]}<b>${c.freq[k] == null ? "not detected" : pct(10 ** c.freq[k], 2)}</b></div>`).join(""), ev);
      })
      .on("mouseleave", () => { paths.style("opacity", 1); TIL.hideTip(); });
  }
  $("#patient").addEventListener("change", drawSlope);

  const all = () => { drawUmap(); drawStrip(); drawPooled(); drawSlope(); };
  let t; window.addEventListener("resize", () => { clearTimeout(t); t = setTimeout(all, 150); });
  window.addEventListener("themechange", all);
  const first = index.datasets.find((d) => !d.synthetic) || index.datasets[0];
  $("#dataset").value = first.id;
  select(first.id);
})();
