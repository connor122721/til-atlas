// Trial tracker page. Data: data/trials.json (written by scraper/fetch_trials.py).
(async function () {
  const { fmt, esc, css } = TIL;
  const $ = (s) => document.querySelector(s);

  const { data } = await TIL.load("data/trials.json");
  const log = (await TIL.load("data/changelog.json")).data || [];
  if (!data) { TIL.emptyState($("#main"), "The first weekly pull hasn't run yet",
    `This page fills itself from the ClinicalTrials.gov API every Monday. To fill it now, open the repository's <b>Actions</b> tab,
     choose <b>Scrape ClinicalTrials.gov</b>, and click <b>Run workflow</b>; the site redeploys about a minute later.`); return; }

  const trials = data.trials;
  $("#meta").textContent = `${fmt(trials.length)} trials · pulled ${data.meta.generated} from ${data.meta.source}`;

  const STATUS_LABEL = (s) => (s || "UNKNOWN").toLowerCase().replace(/_/g, " ").replace(/^./, (c) => c.toUpperCase());
  const PHASE_LABEL = (p) => ({ EARLY_PHASE1: "Early Ph 1", PHASE1: "Ph 1", PHASE2: "Ph 2", PHASE3: "Ph 3", PHASE4: "Ph 4", NA: "N/A" }[p] || p);
  const phaseStr = (t) => t.phases.length ? t.phases.map(PHASE_LABEL).join("/") : "—";
  const year = (d) => (d ? +d.slice(0, 4) : null);

  // ---------- build filter options from the data itself
  const uniq = (arr) => [...new Set(arr)].filter(Boolean);
  const opt = (sel, values, label = (v) => v) => values.forEach((v) => sel.insertAdjacentHTML("beforeend", `<option value="${esc(v)}">${esc(label(v))}</option>`));
  opt($("#status"), uniq(trials.map((t) => t.status)).sort(), STATUS_LABEL);
  opt($("#phase"), uniq(trials.map(phaseStr)).sort());
  opt($("#topic"), data.meta.topics.map((t) => t.id), (id) => data.meta.topics.find((t) => t.id === id).label);
  const designTags = uniq(trials.flatMap((t) => [...(t.tags.combination || []), ...(t.tags.engineering || [])])).sort();
  opt($("#combo"), designTags);
  const indCounts = d3.rollups(trials.flatMap((t) => t.tags.indication), (v) => v.length, (d) => d).sort((a, b) => b[1] - a[1]);
  const selectedInd = new Set();
  $("#indications").innerHTML = indCounts.map(([k, n]) =>
    `<button type="button" class="chip" aria-pressed="false" data-ind="${esc(k)}">${esc(k)} · ${n}</button>`).join("");
  $("#indications").addEventListener("click", (e) => {
    const b = e.target.closest("button[data-ind]"); if (!b) return;
    const on = b.getAttribute("aria-pressed") !== "true";
    b.setAttribute("aria-pressed", on); on ? selectedInd.add(b.dataset.ind) : selectedInd.delete(b.dataset.ind);
    update();
  });

  // ---------- filtering
  const searchText = (t) => [t.id, t.title, t.acronym, t.sponsor, ...t.conditions, ...t.interventions].join(" ").toLowerCase();
  trials.forEach((t) => (t._s = searchText(t)));
  function filtered() {
    const q = $("#q").value.trim().toLowerCase(), st = $("#status").value, ph = $("#phase").value,
      tp = $("#topic").value, cb = $("#combo").value;
    return trials.filter((t) =>
      (!q || q.split(/\s+/).every((w) => t._s.includes(w))) &&
      (!st || t.status === st) && (!ph || phaseStr(t) === ph) && (!tp || t.topics.includes(tp)) &&
      (!cb || [...(t.tags.combination || []), ...(t.tags.engineering || [])].includes(cb)) &&
      (!selectedInd.size || t.tags.indication.some((i) => selectedInd.has(i))));
  }

  // ---------- stat tiles
  function stats(rows) {
    const recruiting = rows.filter((t) => t.status === "RECRUITING").length;
    const countries = uniq(rows.flatMap((t) => t.locations.map((l) => l.country))).length;
    const enroll = rows.map((t) => t.enrollment).filter((n) => n != null);
    const thisYear = new Date().getFullYear();
    const recent = rows.filter((t) => year(t.firstPosted) >= thisYear - 1).length;
    $("#stats").innerHTML = [
      ["Matching trials", fmt(rows.length), `of ${fmt(trials.length)} tracked`],
      ["Recruiting now", fmt(recruiting), rows.length ? `${Math.round(100 * recruiting / rows.length)}% of matches` : ""],
      ["Countries with a site", fmt(countries), ""],
      [`Registered since ${thisYear - 1}`, fmt(recent), enroll.length ? `median enrollment ${fmt(d3.median(enroll))}` : ""],
    ].map(([l, v, n]) => `<div class="card stat"><div class="label">${l}</div><div class="value">${v}</div><div class="note">${n}</div></div>`).join("");
  }

  // ---------- charts
  let yearData = [], indData = [];
  const drawYear = TILCharts.columns($("#byYear"), yearData, { unit: "trials" });
  const drawInd = TILCharts.bars($("#byIndication"), indData, { unit: "trials" });
  function charts(rows) {
    const ys = rows.map((t) => year(t.firstPosted)).filter(Boolean);
    const lo = d3.min(ys), hi = d3.max(ys);
    yearData.length = 0;
    if (ys.length) for (let y = lo; y <= hi; y++) {
      const n = ys.filter((v) => v === y).length;
      yearData.push({ label: String(y), value: n, tip: `${n} trial${n === 1 ? "" : "s"} first posted` });
    }
    indData.length = 0;
    d3.rollups(rows.flatMap((t) => t.tags.indication), (v) => v.length, (d) => d)
      .sort((a, b) => b[1] - a[1]).forEach(([label, value]) => indData.push({ label, value }));
    drawYear(); drawInd();
  }

  // ---------- map
  const map = L.map("map", { scrollWheelZoom: false, worldCopyJump: true }).setView([30, 0], 2);
  let tiles;
  function setTiles() {
    const dark = getComputedStyle(document.documentElement).colorScheme === "dark";
    if (tiles) map.removeLayer(tiles);
    tiles = L.tileLayer(`https://{s}.basemaps.cartocdn.com/${dark ? "dark_all" : "light_all"}/{z}/{x}/{y}{r}.png`, {
      attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> &copy; <a href="https://carto.com/">CARTO</a>',
      subdomains: "abcd", maxZoom: 12 }).addTo(map);
  }
  setTiles(); window.addEventListener("themechange", setTiles);
  const layer = L.layerGroup().addTo(map);
  function drawMap(rows) {
    layer.clearLayers();
    const byCity = new Map();
    rows.forEach((t) => {
      const seen = new Set();
      t.locations.forEach((l) => {
        if (l.lat == null) return;
        const k = `${l.city}|${l.country}`; if (seen.has(k)) return; seen.add(k);
        if (!byCity.has(k)) byCity.set(k, { city: l.city, country: l.country, lat: l.lat, lon: l.lon, trials: [] });
        byCity.get(k).trials.push(t);
      });
    });
    const max = d3.max([...byCity.values()], (c) => c.trials.length) || 1;
    byCity.forEach((c) => {
      L.circleMarker([c.lat, c.lon], { radius: 4 + 12 * Math.sqrt(c.trials.length / max), color: css("--surface-1"), weight: 2,
        fillColor: css("--s1"), fillOpacity: 0.75 })
        .bindTooltip(`<b>${esc(c.city)}, ${esc(c.country)}</b><br>${c.trials.length} trial${c.trials.length > 1 ? "s" : ""}<br>` +
          c.trials.slice(0, 5).map((t) => `<span class="mono">${t.id}</span>`).join("<br>") + (c.trials.length > 5 ? "<br>…" : ""))
        .addTo(layer);
    });
  }

  // ---------- changelog
  function drawLog() {
    const byId = new Map(trials.map((t) => [t.id, t]));
    const e = log[0];
    if (!e) { $("#changelog").innerHTML = `<p class="muted">No pulls recorded yet.</p>`; return; }
    $("#logDate").textContent = `Pull of ${e.date}` + (log[1] ? ` vs ${log[1].date}` : "");
    if (e.initial) { $("#changelog").innerHTML = `<p class="muted">First pull — changes appear from next week on.</p>`; return; }
    const link = (id) => `<a class="mono" href="https://clinicaltrials.gov/study/${id}" target="_blank" rel="noopener">${id}</a>`;
    const parts = [];
    parts.push(`<p><b>${e.new.length}</b> new · <b>${e.statusChanges.length}</b> status changes · <b>${e.removed.length}</b> dropped out of the search</p>`);
    if (e.new.length) parts.push(`<h3 style="margin-top:12px">New</h3>` + e.new.slice(0, 8).map((id) =>
      `<p style="margin:0 0 6px">${link(id)} ${esc(byId.get(id)?.title || "")}</p>`).join(""));
    if (e.statusChanges.length) parts.push(`<h3 style="margin-top:12px">Status changes</h3>` + e.statusChanges.slice(0, 8).map((c) =>
      `<p style="margin:0 0 6px">${link(c.id)}<br><span class="status ${c.from}">${STATUS_LABEL(c.from)}</span> → <span class="status ${c.to}">${STATUS_LABEL(c.to)}</span></p>`).join(""));
    $("#changelog").innerHTML = parts.join("");
  }

  // ---------- table
  let sortKey = "start", sortDir = "descending", shown = 50, current = [];
  const keyFn = {
    id: (t) => t.id, title: (t) => t.title, status: (t) => t.status, phase: phaseStr, indication: (t) => t.tags.indication.join(),
    sponsor: (t) => t.sponsor || "", start: (t) => t.start || "", enrollment: (t) => t.enrollment ?? -1,
  };
  function drawTable() {
    const f = keyFn[sortKey], dir = sortDir === "ascending" ? 1 : -1;
    const rows = [...current].sort((a, b) => (f(a) > f(b) ? 1 : f(a) < f(b) ? -1 : 0) * dir);
    $("#table tbody").innerHTML = rows.slice(0, shown).map((t) => `<tr>
      <td><a class="mono" href="https://clinicaltrials.gov/study/${t.id}" target="_blank" rel="noopener">${t.id}</a></td>
      <td>${esc(t.title)}${t.hasResults ? ' <span class="chip">results posted</span>' : ""}</td>
      <td><span class="status ${t.status}">${STATUS_LABEL(t.status)}</span></td>
      <td>${phaseStr(t)}</td><td>${esc(t.tags.indication.join(", "))}</td>
      <td>${esc(t.sponsor || "")}</td><td class="mono">${t.start || "—"}</td><td class="num">${fmt(t.enrollment)}</td></tr>`).join("");
    $("#more").hidden = rows.length <= shown;
    $("#tableTitle").textContent = `Trials (${fmt(rows.length)})`;
  }
  TIL.sortable($("#table"), (k, d) => { sortKey = k; sortDir = d; drawTable(); });
  $("#table th[data-key=start]").setAttribute("aria-sort", "descending");
  $("#more").addEventListener("click", () => { shown += 100; drawTable(); });

  // ---------- CSV export of exactly what's filtered
  $("#csv").addEventListener("click", () => {
    const cols = ["id", "title", "status", "phase", "indication", "sponsor", "sponsorClass", "start", "firstPosted", "enrollment", "interventions", "countries"];
    const q = (v) => `"${String(v ?? "").replace(/"/g, '""')}"`;
    const lines = [cols.join(",")].concat(current.map((t) => [t.id, t.title, t.status, phaseStr(t), t.tags.indication.join("; "), t.sponsor,
      t.sponsorClass, t.start, t.firstPosted, t.enrollment, t.interventions.join("; "), uniq(t.locations.map((l) => l.country)).join("; ")].map(q).join(",")));
    const a = Object.assign(document.createElement("a"), { href: URL.createObjectURL(new Blob([lines.join("\n")], { type: "text/csv" })),
      download: `til-trials-${data.meta.generated}.csv` });
    a.click(); URL.revokeObjectURL(a.href);
  });

  function update() { current = filtered(); shown = 50; stats(current); charts(current); drawMap(current); drawTable(); }
  ["#q", "#status", "#phase", "#topic", "#combo"].forEach((s) => $(s).addEventListener(s === "#q" ? "input" : "change", update));
  $("#reset").addEventListener("click", () => {
    ["#q", "#status", "#phase", "#topic", "#combo"].forEach((s) => ($(s).value = ""));
    selectedInd.clear(); document.querySelectorAll("#indications .chip").forEach((b) => b.setAttribute("aria-pressed", "false"));
    update();
  });
  drawLog(); update();
})();
