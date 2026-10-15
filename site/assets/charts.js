// Small reusable D3 charts. Thin marks, 4px rounded data-end anchored on the baseline,
// recessive axes, hover tooltip on every mark. Redraw on resize and theme change.
(function () {
  const { css, fmt, esc } = window.TIL;

  /** Path for a bar whose data-end corners are rounded (r) and baseline corners square. */
  function barPath(x, y, w, h, r, horizontal) {
    if (w <= 0 || h <= 0) return "";
    if (horizontal) { r = Math.min(r, h / 2, w);
      return `M${x},${y}H${x + w - r}Q${x + w},${y} ${x + w},${y + r}V${y + h - r}Q${x + w},${y + h} ${x + w - r},${y + h}H${x}Z`; }
    r = Math.min(r, w / 2, h);
    return `M${x},${y + h}V${y + r}Q${x},${y} ${x + r},${y}H${x + w - r}Q${x + w},${y} ${x + w},${y + r}V${y + h}Z`;
  }

  function onRedraw(el, draw) {
    let w = 0;
    const ro = new ResizeObserver(() => { if (el.clientWidth !== w) { w = el.clientWidth; draw(); } });
    ro.observe(el);
    window.addEventListener("themechange", draw);
    return draw;
  }

  /** Vertical bars. data: [{label, value, tip?}] */
  function columns(el, data, { height = 220, color = "--s1", unit = "trials" } = {}) {
    const draw = () => {
      el.innerHTML = "";
      if (!data.length) { el.innerHTML = `<p class="muted small">No trials match these filters.</p>`; return; }
      const W = el.clientWidth || 500, H = height, m = { t: 12, r: 6, b: 26, l: 34 };
      const x = d3.scaleBand().domain(data.map((d) => d.label)).range([m.l, W - m.r]).paddingInner(0.18).paddingOuter(0.05);
      const y = d3.scaleLinear().domain([0, d3.max(data, (d) => d.value) || 1]).nice().range([H - m.b, m.t]);
      const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img")
        .attr("aria-label", `Bar chart of ${unit}`);
      svg.append("g").attr("class", "gridline").attr("transform", `translate(${m.l},0)`)
        .call(d3.axisLeft(y).ticks(4).tickSize(-(W - m.l - m.r)).tickFormat(""));
      svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`)
        .call(d3.axisLeft(y).ticks(4).tickSize(0).tickPadding(6)).call((g) => g.select(".domain").remove());
      const every = Math.ceil(data.length / Math.max(2, Math.floor((W - m.l) / 44)));
      svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`)
        .call(d3.axisBottom(x).tickSize(0).tickPadding(7).tickValues(x.domain().filter((_, i) => i % every === 0 || i === data.length - 1)));
      const bw = Math.min(x.bandwidth(), 28);
      const g = svg.append("g");
      g.selectAll("path").data(data).join("path")
        .attr("d", (d) => barPath(x(d.label) + (x.bandwidth() - bw) / 2, y(d.value), bw, y(0) - y(d.value), 4))
        .attr("fill", css(color));
      // invisible full-height hit targets, larger than the mark
      g.selectAll("rect.hit").data(data).join("rect").attr("class", "hit")
        .attr("x", (d) => x(d.label)).attr("width", x.bandwidth()).attr("y", m.t).attr("height", H - m.b - m.t)
        .attr("fill", "transparent")
        .on("mousemove", (ev, d) => TIL.showTip(`<div class="tt-title">${esc(d.label)}</div>${d.tip || `${fmt(d.value)} ${unit}`}`, ev))
        .on("mouseleave", TIL.hideTip);
      svg.append("line").attr("x1", m.l).attr("x2", W - m.r).attr("y1", y(0)).attr("y2", y(0)).attr("stroke", css("--axis"));
    };
    return onRedraw(el, draw);
  }

  /** Horizontal bars with the value printed at the bar end. data: [{label, value}] */
  function bars(el, data, { color = "--s1", unit = "trials", rowH = 24 } = {}) {
    const draw = () => {
      el.innerHTML = "";
      if (!data.length) { el.innerHTML = `<p class="muted small">No trials match these filters.</p>`; return; }
      const W = el.clientWidth || 500, labelW = Math.min(130, W * 0.36), m = { t: 4, r: 36, b: 4, l: labelW };
      const H = data.length * rowH + m.t + m.b;
      const y = d3.scaleBand().domain(data.map((d) => d.label)).range([m.t, H - m.b]).paddingInner(0.25);
      const x = d3.scaleLinear().domain([0, d3.max(data, (d) => d.value) || 1]).range([m.l, W - m.r]);
      const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", `Bar chart of ${unit}`);
      const row = svg.selectAll("g").data(data).join("g");
      row.append("text").attr("x", m.l - 8).attr("y", (d) => y(d.label) + y.bandwidth() / 2).attr("dy", "0.35em")
        .attr("text-anchor", "end").attr("fill", css("--text-secondary")).attr("font-size", 12.5).text((d) => d.label);
      row.append("path").attr("d", (d) => barPath(m.l, y(d.label), x(d.value) - m.l, y.bandwidth(), 4, true)).attr("fill", css(color));
      row.append("text").attr("x", (d) => x(d.value) + 6).attr("y", (d) => y(d.label) + y.bandwidth() / 2).attr("dy", "0.35em")
        .attr("fill", css("--text-secondary")).attr("font-size", 12).style("font-variant-numeric", "tabular-nums").text((d) => fmt(d.value));
      row.append("rect").attr("x", 0).attr("width", W).attr("y", (d) => y(d.label) - 2).attr("height", y.bandwidth() + 4).attr("fill", "transparent")
        .on("mousemove", (ev, d) => TIL.showTip(`<div class="tt-title">${esc(d.label)}</div>${fmt(d.value)} ${unit}`, ev))
        .on("mouseleave", TIL.hideTip);
      svg.append("line").attr("x1", m.l).attr("x2", m.l).attr("y1", m.t).attr("y2", H - m.b).attr("stroke", css("--axis"));
    };
    return onRedraw(el, draw);
  }

  window.TILCharts = { barPath, columns, bars, onRedraw };
})();
