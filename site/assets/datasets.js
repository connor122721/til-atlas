// Dataset catalog (data/catalog.json, curated) + GEO watch (data/geo_watch.json, weekly).
(async function () {
  const { esc, fmt } = TIL;
  const $ = (s) => document.querySelector(s);
  const cat = (await TIL.load("data/catalog.json")).data;
  const watch = await TIL.load("data/geo_watch.json");
  if (!cat) { TIL.banner($("#main"), "Catalog missing.", "Expected <code>data/catalog.json</code>."); return; }

  const ACCESS = { open: ["Open", "good"], check: ["Check release", "warning"], unknown: ["Not confirmed", "muted"], none: ["No public data", "critical"] };
  const USE = { lab: ["Single-cell", "lab.html"], bulk: ["Bulk RNA", "bulk.html"], spatial: ["Spatial", "spatial.html"], sim: ["Clone simulator", "simulator.html"] };
  const modKey = (m) => (/spatial|xenium|cosmx|visium|imaging|geomx|multiplex/i.test(m) ? "Spatial / imaging"
    : /tcr|bcr/i.test(m) ? "TCR / BCR repertoire" : /scrna|smart-seq/i.test(m) ? "Single-cell RNA" : /bulk/i.test(m) ? "Bulk RNA" : "Other");
  cat.datasets.forEach((d) => (d._mods = [...new Set(d.modalities.map(modKey))]));

  $("#catMeta").textContent = `${cat.datasets.length} datasets · ${cat.datasets.filter((d) => d.onSite).length} on the site now · ${cat.datasets.filter((d) => d.access === "open").length} open · curated ${cat.updated}`;
  [...new Set(cat.datasets.map((d) => d.therapy))].sort().forEach((t) => $("#fTherapy").insertAdjacentHTML("beforeend", `<option>${esc(t)}</option>`));
  const mods = [...new Set(cat.datasets.flatMap((d) => d._mods))].sort();
  const picked = new Set();
  $("#fMod").innerHTML = mods.map((m) => `<button type="button" class="chip" aria-pressed="false" data-m="${esc(m)}">${esc(m)}</button>`).join("");
  $("#fMod").addEventListener("click", (e) => {
    const b = e.target.closest("[data-m]"); if (!b) return;
    const on = b.getAttribute("aria-pressed") !== "true"; b.setAttribute("aria-pressed", on);
    on ? picked.add(b.dataset.m) : picked.delete(b.dataset.m); draw();
  });
  ["#fTherapy", "#fAccess", "#fUse"].forEach((s) => $(s).addEventListener("change", draw));

  const pill = (a, on) => (on ? `<span class="pill on">On the site</span> ` : "") + `<span class="pill ${ACCESS[a][1]}">${ACCESS[a][0]}</span>`;
  function draw() {
    const th = $("#fTherapy").value, ac = $("#fAccess").value, us = $("#fUse").value;
    const rows = cat.datasets.filter((d) => (!th || d.therapy === th) && (!ac || d.access === ac) && (!us || d.use.includes(us)) &&
      (!picked.size || d._mods.some((m) => picked.has(m))));
    $("#cards").innerHTML = rows.map((d) => `<article class="card dataset">
        <div class="ds-head"><h3>${esc(d.short)}</h3><div class="pills">${pill(d.access, d.onSite)}</div></div>
        <p class="small muted">${esc(d.paper)}</p>
        <div class="chips"><span class="chip">${esc(d.cancer)}</span><span class="chip">${esc(d.therapy)}</span>
          ${d.trial ? `<span class="chip mono">${esc(d.trial)}</span>` : ""}${d.modalities.map((m) => `<span class="chip">${esc(m)}</span>`).join("")}</div>
        <p class="small" style="margin-top:10px">${esc(d.samples)}</p>
        <ul class="small acc">${d.accessions.map((a) => `<li>${a.url ? `<a href="${esc(a.url)}" target="_blank" rel="noopener">${esc(a.label)}</a>` : esc(a.label)}</li>`).join("") || "<li class='muted'>No public accession found</li>"}</ul>
        <p class="small muted">${esc(d.status)}</p>
        <div class="ds-foot small">${d.use.map((u) => `<a href="${USE[u][1]}">${USE[u][0]}</a>`).join(" · ")}
          ${d.pipeline ? `<span class="muted"> · pipeline <code>${esc(d.pipeline)}</code></span>` : `<span class="muted"> · no pipeline yet</span>`}</div>
      </article>`).join("") || `<p class="muted">No datasets match.</p>`;
  }
  draw();

  // ---------- GEO watch feed
  const w = watch.data;
  if (!w) { $("#geoMeta").textContent = "The weekly GEO search runs with the trial scraper; new series appear here after its first run."; return; }
  $("#geoMeta").innerHTML = `Human GEO series from the last ${w.meta.lookbackDays} days matching ${w.meta.queries.map((q) => `<i>${esc(q.label)}</i>`).join(", ")}. ` +
    `Last run ${w.meta.generated}: ${w.meta.newThisRun} new. Already-cataloged series are marked.`;
  const items = w.items;
  $("#geo").insertAdjacentHTML("beforeend", `<div class="table-wrap"><table><thead><tr>
      <th>Series</th><th>Title</th><th>Data</th><th class="num">Samples</th><th>Published</th><th>First seen</th></tr></thead><tbody>${
    items.map((i) => `<tr><td><a class="mono" href="${esc(i.url)}" target="_blank" rel="noopener">${esc(i.accession)}</a>${i.inCatalog ? '<br><span class="pill good">in catalog</span>' : ""}</td>
      <td>${esc(i.title)}<div class="small muted">${esc(i.summary.slice(0, 220))}${i.summary.length > 220 ? "…" : ""}</div></td>
      <td><div class="chips">${[...i.modalities, ...i.therapies].map((m) => `<span class="chip">${esc(m)}</span>`).join("")}</div></td>
      <td class="num">${fmt(i.samples)}</td><td class="mono">${esc(i.published)}</td><td class="mono">${esc(i.firstSeen)}</td></tr>`).join("")
  }</tbody></table></div>`);
})();
