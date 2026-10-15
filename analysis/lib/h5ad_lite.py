"""Minimal readers for .h5ad files with h5py only (no anndata/scanpy needed, low memory).

Handles both old-style (obs/__categories) and new-style (obs/<col>/{categories,codes}) categoricals,
and CSR/CSC sparse X. Category code -1 means missing and is returned as None.
"""
from __future__ import annotations

import h5py
import numpy as np
import scipy.sparse as sp


def _dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a], dtype=object)


def obs_col(f: h5py.File, name: str, group: str = "obs") -> np.ndarray:
    g = f[group]
    item = g[name]
    if isinstance(item, h5py.Group):                      # new-style categorical
        cats, codes = _dec(item["categories"][:]), item["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], None)
    v = item[:]
    if "__categories" in g and name in g["__categories"]:  # old-style categorical
        cats = _dec(g["__categories"][name][:])
        return np.where(v >= 0, cats[np.clip(v, 0, None)], None)
    return _dec(v) if v.dtype.kind in "SO" else v


def var_names(f: h5py.File) -> np.ndarray:
    v = f["var"]
    key = v.attrs.get("_index", "_index")
    return _dec(v[key][:])


def obs_names(f: h5py.File) -> np.ndarray:
    o = f["obs"]
    key = o.attrs.get("_index", "_index")
    return _dec(o[key][:])


def read_X(f: h5py.File, key: str = "X") -> sp.csr_matrix:
    X = f[key]
    if isinstance(X, h5py.Dataset):
        return sp.csr_matrix(X[:])
    enc = X.attrs.get("encoding-type", "csr_matrix")
    enc = enc.decode() if isinstance(enc, bytes) else enc
    shape = tuple(X.attrs.get("shape", X.attrs.get("h5sparse_shape")))
    M = (sp.csr_matrix if "csr" in enc else sp.csc_matrix)((X["data"][:], X["indices"][:], X["indptr"][:]), shape=shape)
    return M.tocsr()
