(async function () {
  const { css, fmt } = TIL;
  const $ = (s) => document.querySelector(s);

  // hero numbers from the latest trial pull (demo fixture until the first Action run)
  const { data, path } = await TIL.load("data/trials.json", "data/demo/trials.json");
  if (data) {
    const t = data.trials, demo = path.includes("/demo/");
    const countries = new Set(t.flatMap((x) => x.locations.map((l) => l.country))).size;
    $("#heroStats").innerHTML = [
      ["Trials tracked", fmt(t.length)], ["Recruiting", fmt(t.filter((x) => x.status === "RECRUITING").length)], ["Countries", fmt(countries)],
    ].map(([l, v]) => `<div class="card stat"><div class="label">${l}</div><div class="value">${v}</div></div>`).join("") +
      `<p class="small muted" style="grid-column:1/-1;margin:0">${demo ? "Demo fixture — real numbers appear after the first weekly scrape." : `As of ${data.meta.generated}.`}</p>`;
  }

  // hero art: Muller plot of one responder-like simulated patient
  function art() {
    const el = $("#heroArt"); el.innerHTML = "";
    const sim = SimModel.simulate({ seed: 11, stemFrac: 0.35, reactiveFrac: 0.45, avidity: 1.2, nClones: 60, days: 180 });
    const tl = sim.timeline, cl = sim.clones;
    const rank = { high: 0, low: 1, bystander: 2 };
    const keys = cl.map((_, i) => i).sort((a, b) => rank[cl[a].cls] - rank[cl[b].cls] || cl[b].infused - cl[a].infused);
    const rows = tl.t.filter((_, i) => i % 3 === 0).map((t, j) => { const o = { t }; keys.forEach((k) => (o[k] = Math.sqrt(tl.clones[k][j * 3]))); return o; }); // sqrt keeps the long tail visible
    const st = d3.stack().keys(keys).offset(d3.stackOffsetSilhouette).order(d3.stackOrderInsideOut)(rows);
    const W = el.clientWidth, H = el.clientHeight || W / 1.25;
    const x = d3.scaleLinear().domain([0, rows[rows.length - 1].t]).range([0, W]);
    const y = d3.scaleLinear().domain([d3.min(st, (s) => d3.min(s, (d) => d[0])), d3.max(st, (s) => d3.max(s, (d) => d[1]))]).range([H - 4, 4]);
    const color = { high: "--s1", low: "--s3", bystander: "--s2" };
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
    svg.selectAll("path").data(st).join("path")
      .attr("d", d3.area().curve(d3.curveBasis).x((d) => x(d.data.t)).y0((d) => y(d[0])).y1((d) => y(d[1])))
      .attr("fill", (s) => css(color[cl[s.key].cls])).attr("stroke", css("--surface-1")).attr("stroke-width", 0.5);
  }

  // pipeline diagram
  function flow() {
    const el = $("#flow"); el.innerHTML = "";
    const steps = [["Weekly cron", "GitHub Action, Mondays"], ["Query API v2", "paginate with nextPageToken"], ["Normalize + tag", "regex rules in tags.json"],
      ["Diff vs last week", "new · status changes · removed"], ["Commit JSON", "site/data/*.json"], ["Pages deploy", "static site rebuilds"]];
    const W = el.clientWidth, narrow = W < 640, cols = narrow ? 2 : 6, bw = (W - (cols - 1) * 22) / cols, bh = 62, rowsN = Math.ceil(steps.length / cols);
    const H = rowsN * bh + (rowsN - 1) * 24 + 4;
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img").attr("aria-label", "Pipeline: cron, API, normalize, diff, commit, deploy");
    steps.forEach(([a, b], i) => {
      const cx = (i % cols) * (bw + 22), cy = Math.floor(i / cols) * (bh + 24) + 2;
      svg.append("rect").attr("x", cx).attr("y", cy).attr("width", bw).attr("height", bh).attr("rx", 8)
        .attr("fill", css("--surface-2")).attr("stroke", css("--border"));
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
