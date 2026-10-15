// Bulk RNA tab: forest plot of AUC by cohort, per-cohort strip plot, Kaplan–Meier, AUC heatmap.
// Data: data/bulk/index.json, summary.json, <cohort>.json (analysis/bulk/fetch_iatlas.py).
(async function () {
  const { css, fmt, pct, esc } = TIL, { rankSumP, pTxt, strip } = TILStats;
  const $ = (s) => document.querySelector(s);
  const idx = (await TIL.load("data/bulk/index.json")).data;
  const sum = (await TIL.load("data/bulk/summary.json")).data;
  if (!idx) { TIL.banner($("#main"), "No bulk data yet.", "Run <code>python analysis/bulk/fetch_iatlas.py</code>."); return; }
  const S = sum.stats, cohorts = idx.cohorts, cache = {};
  const totalN = d3.sum(cohorts, (c) => c.n);
  $("#meta").textContent = `${cohorts.length} cohorts · ${fmt(totalN)} patients · ${[...new Set(cohorts.map((c) => c.cancer))].join(", ")}`;

  const computed = Object.keys(idx.signatures).filter((k) => !k.startsWith("ia_"));
  const iatlas = Object.keys(idx.signatures).filter((k) => k.startsWith("ia_"));
  $("#sig").innerHTML = `<optgroup label="Computed here">${computed.map((k) => `<option value="${k}">${esc(idx.signatures[k])}</option>`).join("")}</optgroup>` +
    `<optgroup label="Precomputed by iAtlas">${iatlas.map((k) => `<option value="${k}">${esc(idx.signatures[k])}</option>`).join("")}</optgroup>`;
  $("#ep").innerHTML = Object.entries(idx.endpoints).map(([k, v]) => `<option value="${k}">${esc(v)}</option>`).join("");
  $("#gene").innerHTML = `<option value="">Signature score</option>` + idx.panel.map((g) => `<option value="${g}">Gene: ${g}</option>`).join("");
  let cur = cohorts[0].id;
  const load = async (id) => (cache[id] ||= (await TIL.load(`data/bulk/${id}.json`)).data);
  const sig = () => $("#sig").value, ep = () => $("#ep").value;

  // ---------------- forest plot
  function drawForest() {
    const el = $("#forest"); el.innerHTML = "";
    const k = sig(), e = ep();
    $("#forestTitle").textContent = `${idx.signatures[k]} vs ${e === "response" ? "objective response" : "clinical benefit"}`;
    $("#forestSub").textContent = `${idx.signatureSources[k]}. AUC with 95% bootstrap interval; dot size = patients. Click a cohort.`;
    const rows = cohorts.map((c) => ({ c, s: S[c.id]?.[e]?.[k] })).filter((r) => r.s);
    const W = el.clientWidth || 640, rowH = 34, lab = Math.min(230, W * 0.42), m = { t: 22, r: 16, b: 30, l: lab };
    const H = m.t + rows.length * rowH + m.b;
    const x = d3.scaleLinear().domain([0.15, 1]).range([m.l, W - m.r]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Forest plot of AUC by cohort");
    svg.append("g").attr("class", "gridline").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([0.25, 0.5, 0.75, 1]).tickSize(-(H - m.t - m.b)).tickFormat(""));
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([0.25, 0.5, 0.75, 1]).tickFormat(d3.format(".2~f")).tickSize(0).tickPadding(8)).call((g) => g.select(".domain").remove());
    svg.append("line").attr("x1", x(0.5)).attr("x2", x(0.5)).attr("y1", m.t - 8).attr("y2", H - m.b).attr("stroke", css("--text-muted")).attr("stroke-dasharray", "4 3");
    svg.append("text").attr("x", x(0.5) - 6).attr("y", m.t - 10).attr("text-anchor", "end").attr("font-size", 11).attr("fill", css("--text-muted")).text("← non-responders higher");
    svg.append("text").attr("x", x(0.5) + 6).attr("y", m.t - 10).attr("font-size", 11).attr("fill", css("--text-muted")).text("responders higher →");
    const r = d3.scaleSqrt().domain([0, d3.max(rows, (d) => d.s.nR + d.s.nNR)]).range([3, 9]);
    const g = svg.selectAll("g.row").data(rows).join("g").attr("class", "row").attr("transform", (d, i) => `translate(0,${m.t + i * rowH + rowH / 2})`).style("cursor", "pointer");
    g.append("rect").attr("x", 0).attr("y", -rowH / 2).attr("width", W).attr("height", rowH).attr("fill", (d) => (d.c.id === cur ? css("--surface-2") : "transparent"));
    g.append("text").attr("x", 4).attr("dy", "-0.15em").attr("font-size", 12.5).attr("font-weight", 600).attr("fill", css("--text-primary")).text((d) => d.c.label);
    g.append("text").attr("x", 4).attr("dy", "1.05em").attr("font-size", 11).attr("fill", css("--text-secondary")).text((d) => `${d.c.cancer} · ${d.s.nR} R / ${d.s.nNR} NR`);
    g.append("line").attr("x1", (d) => x(d.s.lo)).attr("x2", (d) => x(d.s.hi)).attr("stroke", css("--s1")).attr("stroke-width", 2);
    g.append("circle").attr("cx", (d) => x(d.s.auc)).attr("r", (d) => r(d.s.nR + d.s.nNR)).attr("fill", css("--s1")).attr("stroke", css("--surface-1")).attr("stroke-width", 2);
    g.on("mousemove", (ev, d) => TIL.showTip(`<div class="tt-title">${esc(d.c.label)} · ${esc(d.c.cancer)}</div>
        <div class="small muted">${esc(d.c.citation)}</div>
        <div class="row">AUC<b>${d.s.auc.toFixed(2)} (${d.s.lo.toFixed(2)}–${d.s.hi.toFixed(2)})</b></div>
        <div class="row">Rank-sum p<b>${pTxt(d.s.p)}</b></div><div class="row">Patients<b>${d.s.nR} R / ${d.s.nNR} NR</b></div>`, ev))
      .on("mouseleave", TIL.hideTip).on("click", (ev, d) => { cur = d.c.id; drawForest(); drawCohort(); });
  }

  // ---------------- per-cohort strip + KM
  async function drawCohort() {
    const D = await load(cur), k = sig(), e = ep(), gene = $("#gene").value;
    const c = D.meta;
    $("#cohortTitle").textContent = `${c.label} · ${c.cancer}`;
    $("#cohortSub").innerHTML = `${esc(c.citation)}${c.pmid ? ` · <a href="https://pubmed.ncbi.nlm.nih.gov/${c.pmid}/" target="_blank" rel="noopener">PubMed</a>` : ""}<br>${esc(c.therapy)}`;
    const val = (s) => (gene ? s.genes[gene] : s.scores[k]);
    const rows = D.samples.filter((s) => s[e] && val(s) != null);
    const groups = ["Responder", "Non-responder"].map((r, i) => ({ key: r, label: r, color: i ? "--s2" : "--s1", rows: rows.filter((s) => s[e] === r) }));
    strip($("#strip"), groups, val, { H: 250, fmtY: d3.format(".1f"), zeroLine: true,
      tip: (s) => `<div class="tt-title">${esc(s.patient)}</div><div class="row">${esc(gene || idx.signatures[k])}<b>${val(s).toFixed(2)}</b></div>
        ${s.recist ? `<div class="row">Best response<b>${esc(s.recist)}</b></div>` : ""}${s.os != null ? `<div class="row">OS<b>${s.os} mo${s.osEvent ? "" : " (alive)"}</b></div>` : ""}` });
    const a = groups[0].rows.map(val), b = groups[1].rows.map(val);
    $("#stripStats").innerHTML = `${esc(gene ? "Gene " + gene + " (z-score)" : idx.signatures[k])}: ${a.length} vs ${b.length} patients, two-sided Wilcoxon <b>p = ${pTxt(rankSumP(a, b))}</b>.`;
    drawKM(D, val, gene || idx.signatures[k]);
  }

  function kmCurve(rows) {
    const t = rows.map((r) => [r.os, r.osEvent]).sort((p, q) => p[0] - q[0]);
    let atRisk = t.length, s = 1; const pts = [[0, 1]];
    for (let i = 0; i < t.length;) {
      const time = t[i][0]; let d = 0, n = 0;
      while (i < t.length && t[i][0] === time) { if (t[i][1]) d++; n++; i++; }
      if (d) { pts.push([time, s]); s *= 1 - d / atRisk; pts.push([time, s]); }
      atRisk -= n;
    }
    pts.push([t[t.length - 1][0], s]);
    return pts;
  }
  function logRankP(g1, g2) { // two-group log-rank, chi-square 1 df
    const all = [...g1.map((r) => [r.os, r.osEvent, 0]), ...g2.map((r) => [r.os, r.osEvent, 1])].sort((a, b) => a[0] - b[0]);
    let n1 = g1.length, n2 = g2.length, O = 0, E = 0, V = 0;
    for (let i = 0; i < all.length;) {
      const t = all[i][0]; let d1 = 0, d = 0, c1 = 0, c = 0;
      while (i < all.length && all[i][0] === t) { const [, ev, g] = all[i]; if (ev) { d++; if (!g) d1++; } c++; if (!g) c1++; i++; }
      const n = n1 + n2;
      if (d && n > 1) { O += d1; E += d * n1 / n; V += d * (n1 / n) * (n2 / n) * (n - d) / (n - 1); }
      n1 -= c1; n2 -= c - c1;
    }
    const chi = V > 0 ? (O - E) ** 2 / V : 0;
    return erfc(Math.sqrt(chi / 2));          // P(chi-square, 1 df > chi)
  }
  function drawKM(D, val, name) {
    const el = $("#km"); el.innerHTML = "";
    const rows = D.samples.filter((s) => s.os != null && s.osEvent != null && val(s) != null);
    if (rows.length < 10) { el.innerHTML = `<p class="muted small">This cohort has no overall survival data.</p>`; $("#kmStats").textContent = ""; $("#kmLegend").innerHTML = ""; return; }
    const med = d3.median(rows, val), hi = rows.filter((s) => val(s) > med), lo = rows.filter((s) => val(s) <= med);
    const W = el.clientWidth || 420, H = 240, m = { t: 10, r: 14, b: 30, l: 40 };
    const x = d3.scaleLinear().domain([0, d3.max(rows, (s) => s.os)]).nice().range([m.l, W - m.r]);
    const y = d3.scaleLinear().domain([0, 1]).range([H - m.b, m.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Kaplan–Meier overall survival");
    svg.append("g").attr("class", "gridline").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(4).tickSize(-(W - m.l - m.r)).tickFormat(""));
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(4, "%").tickSize(0).tickPadding(6)).call((g) => g.select(".domain").remove());
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6).tickSize(0).tickPadding(8)).call((g) => g.select(".domain").remove());
    svg.append("text").attr("x", W - m.r).attr("y", H - 2).attr("text-anchor", "end").attr("font-size", 11).attr("fill", css("--text-muted")).text("months");
    const line = d3.line().x((d) => x(d[0])).y((d) => y(d[1]));
    [[hi, "--s1", "High"], [lo, "--s2", "Low"]].forEach(([g, c]) => svg.append("path").attr("d", line(kmCurve(g))).attr("fill", "none").attr("stroke", css(c)).attr("stroke-width", 2));
    $("#kmLegend").innerHTML = `<span><i style="background:${css("--s1")}"></i>Above median (${hi.length})</span><span><i style="background:${css("--s2")}"></i>At or below median (${lo.length})</span>`;
    const chiP = logRankP(hi, lo);
    $("#kmStats").innerHTML = `${esc(name)}, split at the median. Log-rank <b>p = ${pTxt(chiP)}</b>. Deaths: ${hi.filter((s) => s.osEvent).length} high, ${lo.filter((s) => s.osEvent).length} low.`;
    $("#kmSub").textContent = `${D.meta.label}: ${rows.length} patients with survival data`;
  }
  function erfc(x) { // Abramowitz–Stegun 7.1.26
    const t = 1 / (1 + 0.3275911 * x), y = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    return y * Math.exp(-x * x);
  }

  const fit = (t, n) => (t.length > n ? t.slice(0, Math.max(3, n - 1)) + "…" : t);
  // ---------------- AUC heatmap
  function drawHeat() {
    const el = $("#heat"); el.innerHTML = "";
    const e = ep(), keys = [...computed, ...iatlas], cs = cohorts.filter((c) => S[c.id]?.[e]);
    const W = el.clientWidth || 420, lab = Math.min(150, W * 0.38), top = 92, cw = (W - lab - 4) / cs.length, ch = 18, H = top + keys.length * ch + 4;
    const col = d3.scaleLinear().domain([0.2, 0.35, 0.5, 0.65, 0.8]).range(["--div-p2", "--div-p1", "--div-0", "--div-n1", "--div-n2"].map(css)).clamp(true);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "AUC heatmap");
    cs.forEach((c, j) => svg.append("text").attr("transform", `translate(${lab + j * cw + cw / 2},${top - 6}) rotate(-55)`).attr("font-size", 11).attr("fill", css("--text-secondary")).text(c.label));
    keys.forEach((k, i) => {
      svg.append("text").attr("x", lab - 6).attr("y", top + i * ch + ch / 2).attr("dy", "0.35em").attr("text-anchor", "end").attr("font-size", 11)
        .attr("fill", css(k === sig() ? "--text-primary" : "--text-secondary")).attr("font-weight", k === sig() ? 700 : 400).text(fit(k.startsWith("ia_") ? idx.signatures[k].replace(/ \((.*)\)$/, " · $1") : idx.signatures[k], Math.floor((lab - 8) / 5.6)))
        .append("title").text(idx.signatures[k]);
      cs.forEach((c, j) => {
        const s = S[c.id][e][k];
        svg.append("rect").attr("x", lab + j * cw + 1).attr("y", top + i * ch + 1).attr("width", cw - 2).attr("height", ch - 2).attr("rx", 2)
          .attr("fill", s ? col(s.auc) : css("--surface-2")).style("cursor", "pointer")
          .on("mousemove", (ev) => TIL.showTip(`<div class="tt-title">${esc(idx.signatures[k])}</div><div class="row">${esc(c.label)}<b>${s ? `AUC ${s.auc.toFixed(2)} (${s.lo.toFixed(2)}–${s.hi.toFixed(2)})` : "n/a"}</b></div>`, ev))
          .on("mouseleave", TIL.hideTip)
          .on("click", () => { $("#sig").value = k; cur = c.id; $("#gene").value = ""; all(); });
      });
    });
    const stops = [0.2, 0.35, 0.5, 0.65, 0.8].map(col).join(",");
    $("#hmLegend").innerHTML = `<span>0.2</span><span style="width:120px;height:8px;border-radius:4px;background:linear-gradient(90deg,${stops})"></span><span>0.8 AUC</span>`;
  }

  const all = () => { drawForest(); drawCohort(); drawHeat(); };
  $("#sig").addEventListener("change", () => { $("#gene").value = ""; all(); });
  $("#ep").addEventListener("change", all);
  $("#gene").addEventListener("change", drawCohort);
  let t; window.addEventListener("resize", () => { clearTimeout(t); t = setTimeout(all, 150); });
  window.addEventListener("themechange", all);
  all();
})();
