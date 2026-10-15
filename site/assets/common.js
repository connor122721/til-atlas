// Shared helpers for every page: nav, theme toggle, tooltips, data loading.
(function () {
  const PAGES = [
    ["index.html", "Home"],
    ["trials.html", "Trial tracker"],
    ["simulator.html", "Clone simulator"],
    ["lab.html", "Replication lab"],
    ["spatial.html", "Spatial"],
    ["datasets.html", "Datasets"],
    ["learn.html", "Learn"],
  ];
  const here = location.pathname.split("/").pop() || "index.html";

  const logo = `<svg viewBox="0 0 32 32" aria-hidden="true">
    <circle cx="16" cy="16" r="14" fill="none" stroke="var(--accent)" stroke-width="2"/>
    <circle cx="11" cy="13" r="3.2" fill="var(--s1)"/><circle cx="20" cy="11" r="2.4" fill="var(--s3)"/>
    <circle cx="19" cy="20" r="3.6" fill="var(--s2)"/><circle cx="11" cy="21" r="2" fill="var(--s7)"/></svg>`;

  const header = document.createElement("header");
  header.className = "topbar";
  header.innerHTML = `<div class="inner">
    <a class="brand" href="index.html">${logo}<span>TIL Atlas</span></a>
    <nav class="nav">${PAGES.map(([h, l]) =>
      `<a href="${h}"${h === here ? ' aria-current="page"' : ""}>${l}</a>`).join("")}</nav>
    <button class="theme-btn" type="button" aria-label="Toggle color theme"><svg viewBox="0 0 20 20" width="16" height="16" aria-hidden="true" style="display:block"><circle cx="10" cy="10" r="7.5" fill="none" stroke="currentColor" stroke-width="1.6"/><path d="M10 2.5a7.5 7.5 0 0 1 0 15z" fill="currentColor"/></svg></button></div>`;
  document.body.prepend(header);

  const footer = document.createElement("footer");
  footer.innerHTML = `<div class="inner">TIL Atlas · trial data from the ClinicalTrials.gov API v2, refreshed weekly by a GitHub Action ·
    analyses re-derived from public datasets; see <a href="learn.html">Learn</a> for sources. Not medical advice.</div>`;
  document.body.append(footer);

  // theme: explicit choice beats OS; remembered per browser
  const root = document.documentElement;
  try { const t = localStorage.getItem("til-theme"); if (t) root.dataset.theme = t; } catch (e) {}
  header.querySelector(".theme-btn").addEventListener("click", () => {
    const dark = root.dataset.theme ? root.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
    root.dataset.theme = dark ? "light" : "dark";
    try { localStorage.setItem("til-theme", root.dataset.theme); } catch (e) {}
    window.dispatchEvent(new Event("themechange"));
  });
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => window.dispatchEvent(new Event("themechange")));

  // tooltip singleton
  const tip = document.createElement("div");
  tip.className = "tooltip"; tip.setAttribute("role", "tooltip");
  document.body.append(tip);

  window.TIL = {
    /** Read a CSS custom property, e.g. TIL.css('--s1'). */
    css: (name) => getComputedStyle(root).getPropertyValue(name).trim(),
    fmt: (n, d = 0) => (n == null || isNaN(n) ? "—" : Number(n).toLocaleString(undefined, { maximumFractionDigits: d, minimumFractionDigits: d })),
    pct: (x, d = 0) => (x == null || isNaN(x) ? "—" : (100 * x).toFixed(d) + "%"),
    esc: (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c])),
    showTip(html, ev) {
      tip.innerHTML = html; tip.style.opacity = 1;
      const pad = 14, r = tip.getBoundingClientRect();
      let x = ev.clientX + pad, y = ev.clientY + pad;
      if (x + r.width > innerWidth - 8) x = ev.clientX - r.width - pad;
      if (y + r.height > innerHeight - 8) y = ev.clientY - r.height - pad;
      tip.style.left = Math.max(8, x) + "px"; tip.style.top = Math.max(8, y) + "px";
    },
    hideTip() { tip.style.opacity = 0; },
    /** Fetch JSON; on failure try each fallback path in turn. Returns {data, path}. */
    async load(path, ...fallbacks) {
      for (const p of [path, ...fallbacks]) {
        try {
          const r = await fetch(p, { cache: "no-cache" });
          if (r.ok) return { data: await r.json(), path: p };
        } catch (e) { /* try next */ }
      }
      return { data: null, path: null };
    },
    banner(el, title, body) {
      el.insertAdjacentHTML("afterbegin", `<div class="banner" role="note">
        <svg class="icon" viewBox="0 0 20 20" fill="currentColor" aria-hidden="true"><path d="M10 2 1 18h18L10 2zm0 5 .1 6h-.2L10 7zm0 9.2a1.1 1.1 0 1 1 0-2.2 1.1 1.1 0 0 1 0 2.2z"/></svg>
        <div><b>${title}</b> ${body}</div></div>`);
    },
    /** Make a <table> sortable by clicking headers; numeric when th.num. */
    sortable(table, onSort) {
      table.querySelectorAll("th").forEach((th, i) => th.addEventListener("click", () => {
        const dir = th.getAttribute("aria-sort") === "ascending" ? "descending" : "ascending";
        table.querySelectorAll("th").forEach((o) => o.removeAttribute("aria-sort"));
        th.setAttribute("aria-sort", dir);
        onSort(th.dataset.key || i, dir);
      }));
    },
  };
})();
