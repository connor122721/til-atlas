// Small statistics + a strip plot shared by the lab and spatial pages.
(function () {
  const { css, esc } = TIL;
  const normCdf = (z) => 0.5 * (1 + Math.tanh(Math.sqrt(2 / Math.PI) * (z + 0.044715 * z ** 3)));

  /** Two-sided Wilcoxon rank-sum p: exact (enumerate all relabelings) when n ≤ 22, else normal approximation. */
  function rankSumP(a, b) {
    if (!a.length || !b.length) return NaN;
    const all = [...a, ...b], s = all.map((v, i) => [v, i]).sort((p, q) => p[0] - q[0]), r = Array(all.length);
    for (let i = 0; i < s.length;) { let j = i; while (j + 1 < s.length && s[j + 1][0] === s[i][0]) j++;
      for (let k = i; k <= j; k++) r[s[k][1]] = (i + j) / 2 + 1; i = j + 1; }
    const obs = d3.sum(r.slice(0, a.length)), N = all.length, mean = a.length * (N + 1) / 2;
    if (N > 22) {
      const sd = Math.sqrt(a.length * b.length * (N + 1) / 12);
      return Math.min(1, 2 * (1 - normCdf((Math.abs(obs - mean) - 0.5) / sd)));
    }
    let ext = 0, tot = 0;
    (function rec(st, k, sum) {
      if (k === 0) { tot++; if (Math.abs(sum - mean) >= Math.abs(obs - mean) - 1e-9) ext++; return; }
      for (let i = st; i <= N - k; i++) rec(i + 1, k - 1, sum + r[i]);
    })(0, a.length, 0);
    return ext / tot;
  }
  const pTxt = (p) => (isNaN(p) ? "n/a" : p < 0.001 ? p.toExponential(1) : p.toFixed(3));

  /** Strip plot. groups: [{key,label,color,rows}], value(row) -> number|null, tip(row) -> html. */
  function strip(el, groups, value, { H = 260, fmtY = d3.format(".2~f"), tip = () => "", zeroLine = false } = {}) {
    el.innerHTML = "";
    const vals = groups.flatMap((g) => g.rows.map(value)).filter((v) => v != null && !isNaN(v));
    if (!vals.length) { el.innerHTML = `<p class="muted small">No values for this metric.</p>`; return; }
    const W = el.clientWidth || 300, m = { t: 10, r: 10, b: 28, l: 46 };
    const x = d3.scaleBand().domain(groups.map((g) => g.key)).range([m.l, W - m.r]).padding(0.35);
    let [lo, hi] = d3.extent(vals); if (lo === hi) { lo -= 1; hi += 1; }
    const pad = (hi - lo) * 0.1;
    const y = d3.scaleLinear().domain([Math.min(lo - pad, zeroLine ? 0 : lo - pad), hi + pad]).nice().range([H - m.b, m.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Per-sample values by group");
    svg.append("g").attr("class", "gridline").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickSize(-(W - m.l - m.r)).tickFormat(""));
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickFormat(fmtY).tickSize(0).tickPadding(6)).call((g) => g.select(".domain").remove());
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`)
      .call(d3.axisBottom(x).tickFormat((k) => groups.find((g) => g.key === k).label).tickSize(0).tickPadding(8)).call((g) => g.select(".domain").remove());
    if (zeroLine && y.domain()[0] < 0) svg.append("line").attr("x1", m.l).attr("x2", W - m.r).attr("y1", y(0)).attr("y2", y(0)).attr("stroke", css("--axis"));
    const jr = d3.randomUniform.source(d3.randomLcg(0.41))(-1, 1);
    groups.forEach((g) => {
      const v = g.rows.map(value).filter((u) => u != null && !isNaN(u)); if (!v.length) return;
      const med = d3.median(v);
      svg.append("line").attr("x1", x(g.key)).attr("x2", x(g.key) + x.bandwidth()).attr("y1", y(med)).attr("y2", y(med)).attr("stroke", css("--text-primary")).attr("stroke-width", 2);
      svg.append("g").selectAll("circle").data(g.rows.filter((d) => value(d) != null && !isNaN(value(d)))).join("circle")
        .attr("cx", () => x(g.key) + x.bandwidth() / 2 + jr() * x.bandwidth() * 0.32).attr("cy", (d) => y(value(d))).attr("r", 5)
        .attr("fill", css(g.color)).attr("stroke", css("--surface-1")).attr("stroke-width", 2)
        .on("mousemove", (ev, d) => TIL.showTip(tip(d), ev)).on("mouseleave", TIL.hideTip);
    });
  }
  window.TILStats = { rankSumP, pTxt, strip };
})();
