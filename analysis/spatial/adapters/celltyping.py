"""Marker-module cell typing shared by the Xenium and CosMx adapters (numpy/scipy only)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp

MARKERS = Path(__file__).with_name("markers.json")


def type_cells(X: sp.csr_matrix, genes: list[str], markers: Path = MARKERS) -> tuple[np.ndarray, dict]:
    """X: cells x genes raw counts. Returns (labels, report).

    Per cell: scale to 100 total counts, log1p; per type: mean of its markers; z-score each
    module across cells; label = best module if its z >= threshold, else 'Other'.
    """
    cfg = json.loads(markers.read_text())
    X = sp.csr_matrix(X, dtype=np.float32)
    tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    Xn = sp.diags(100.0 / tot) @ X
    Xn.data = np.log1p(Xn.data)
    gi = {g: i for i, g in enumerate(genes)}
    names, scores, used = [], [], {}
    for t, gs in cfg["types"].items():
        idx = [gi[g] for g in gs if g in gi]
        if len(idx) < 2:
            continue
        s = np.asarray(Xn[:, idx].mean(axis=1)).ravel()
        sd = s.std()
        scores.append((s - s.mean()) / (sd if sd > 0 else 1))
        names.append(t)
        used[t] = [g for g in gs if g in gi]
    if not scores:
        raise SystemExit("none of the marker panels has >= 2 genes in this dataset; edit markers.json")
    S = np.vstack(scores)
    best = S.argmax(axis=0)
    labels = np.where(S.max(axis=0) >= cfg["threshold"], np.array(names)[best], "Other")
    report = {"markersUsed": used, "counts": {t: int((labels == t).sum()) for t in [*names, "Other"]}}
    return labels, report
