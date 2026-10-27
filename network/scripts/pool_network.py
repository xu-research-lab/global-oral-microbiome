"""Pool cohort SparCC matrices and retain reproducible consensus edges."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import numpy as np
import pandas as pd
from scipy import stats
from network_common import MODULE_DIR, RESULT_DIR, GROUPS

ROOT = str(MODULE_DIR)
RESULT = str(RESULT_DIR)
SPARCC_BASE = str(RESULT_DIR / "per_site_sparcc")
RSCRIPT = os.environ.get("RSCRIPT", shutil.which("Rscript") or "Rscript")
NETMOSS_WRAPPER = str(Path(__file__).with_name("netmoss_get_network.R"))
EPS = 1e-8

def bh_fdr(p):
    """Benjamini-Hochberg FDR (NaN preserved)."""
    p = np.asarray(p, dtype=float).ravel()
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return out
    order = np.argsort(pv, kind="mergesort")
    pv = pv[order]
    adj = pv * (n / np.arange(1, n + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    inv = np.empty(n, dtype=np.intp)
    inv[order] = np.arange(n)
    out[ok] = adj[inv]
    return out

def weight_block(r, n):
    """uw: w = (n-1) / (1 - r^2), clipped to avoid divergence as |r| -> 1."""
    return (n - 1) / np.clip(1.0 - r ** 2, EPS, 1.0)

def load_unit_matrix(unit_id, sid, pos):
    fp = os.path.join(SPARCC_BASE, unit_id, f"{sid}.csv")
    if not os.path.exists(fp):
        return None
    df = pd.read_csv(fp, sep=",", index_col=0)
    df.index = [str(x) for x in df.index]
    df.columns = [str(x) for x in df.columns]
    inter = [f for f in df.index if f in pos]
    if len(inter) < 2:
        return None
    df = df.loc[inter, inter]
    idx = np.fromiter((pos[f] for f in inter), dtype=np.int64, count=len(inter))
    return idx, df.to_numpy(dtype=np.float64)

def discover_units(grp, units_all):
    """Units of one group with a SparCC matrix present."""
    units = []
    for _, u in units_all[(units_all.category == grp) & units_all.kept].iterrows():
        fp = os.path.join(SPARCC_BASE, u.unit_id, f"{u.study}.csv")
        if os.path.exists(fp):
            units.append(dict(unit_id=u.unit_id, study=str(u.study),
                              n=int(u.n_samples), path=fp))
        else:
            print(f"[warn] no SparCC matrix for {u.unit_id}/{u.study}, skipped")
    return units

def collect_fid(units):
    """Union genus space across the group's units."""
    fid = set()
    for u in units:
        df = pd.read_csv(u["path"], sep=",", index_col=0, nrows=0)
        fid.update(str(x) for x in df.index)
        fid.update(str(x) for x in df.columns)
    return np.array(sorted(fid))

def build_stats(units, fid):
    F = len(fid)
    pos = {f: k for k, f in enumerate(fid)}

    W1_all = np.zeros((F, F))   # zero-padded weight sum        (consensus)
    VR_all = np.zeros((F, F))   # zero-padded weight*value sum  (consensus)
    W1 = np.zeros((F, F))       # observed weight sum
    VR = np.zeros((F, F))       # observed weight*value sum
    cov = np.zeros((F, F), dtype=np.int32)

    cohorts = []
    for u in units:
        n = u["n"]
        if n < 2:
            continue
        got = load_unit_matrix(u["unit_id"], u["study"], pos)
        if got is None:
            continue
        idx, r = got
        w = weight_block(r, n)

        present = np.zeros(F, dtype=bool)
        present[idx] = True
        obs = np.outer(present, present)
        obsb = obs.astype(np.float64)

        v_full = np.zeros((F, F))
        w_full = np.zeros((F, F))
        v_full[np.ix_(idx, idx)] = r            # uw: pooling space == r space
        w_full[np.ix_(idx, idx)] = w

        w_full[~obs] = n - 1                    # zero-pad missing pairs
        np.fill_diagonal(w_full, 1.0)           # self-correlation weight 1
        np.fill_diagonal(v_full, 0.0)

        W1_all += w_full
        VR_all += w_full * v_full
        W1 += w_full * obsb
        VR += w_full * v_full * obsb
        cov += obs

        cohorts.append((idx, r, w, n, u))
        print(f"  [unit] {u['study']}  n={n}  genera={len(idx)}")

    return dict(W1_all=W1_all, VR_all=VR_all, W1=W1, VR=VR,
                cov=cov), cohorts

def consensus_sign(acc):
    with np.errstate(invalid="ignore", divide="ignore"):
        cons = np.where(acc["W1_all"] > 0, acc["VR_all"] / acc["W1_all"], 0.0)
    cs = np.sign(cons)
    np.fill_diagonal(cs, 1.0)
    return cs

def consensus_pool_py(acc):
    """In-house zero-padded pooling — kept ONLY as the patch source for the
    saturated (non-finite) pairs that NetMoss2 reports as NaN; all results
    come from the official package."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(acc["W1_all"] > 0, acc["VR_all"] / acc["W1_all"], 0.0)

def _manifest_text(units):
    rows = pd.DataFrame([{"study": u["study"],
                          "matrix_path": os.path.abspath(u["path"]),
                          "n": int(u["n"])} for u in units])
    return rows.to_csv(index=False)


def input_fingerprint(units):
    digest = hashlib.sha256()
    for unit in units:
        digest.update(f"{unit['study']}:{unit['n']}:".encode())
        digest.update(Path(unit["path"]).read_bytes())
    return digest.hexdigest()

def netmoss_pooled_matrix(units, fid, rscript, grp):
    """Zero-padded pooled r over fid from the official NetMoss::getNetwork.

    Writes the manifest to its cache dir and reuses the pooled matrix when
    the manifest is unchanged. Returns (cons_matrix F x F, meta dict).
    """
    out_dir = os.path.join(RESULT, "backbone_netmoss", f"{grp}_genus")
    os.makedirs(out_dir, exist_ok=True)
    man = os.path.join(out_dir, "manifest.csv")
    man_text = _manifest_text(units)
    pooled_fp = os.path.join(out_dir, "pooled_union.csv")
    run_meta = {}
    meta_fp = os.path.join(out_dir, "netmoss_run.json")
    if os.path.isfile(meta_fp):
        try:
            with open(meta_fp) as f:
                run_meta = {k: v for k, v in json.load(f).items()}
        except (json.JSONDecodeError, OSError):
            run_meta = {}
    fingerprint = input_fingerprint(units)
    wrapper_fingerprint = hashlib.sha256(Path(NETMOSS_WRAPPER).read_bytes()).hexdigest()
    cached = (os.path.isfile(pooled_fp) and run_meta.get("engine") == "netmoss2"
              and run_meta.get("input_sha256") == fingerprint
              and run_meta.get("wrapper_sha256") == wrapper_fingerprint)
    if not cached:
        with open(man, "w") as f:
            f.write(man_text)
        cmd = [rscript, NETMOSS_WRAPPER, "--manifest", man, "--out", out_dir]
        try:
            r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                               timeout=7200)
        except FileNotFoundError as e:
            raise RuntimeError(f"Rscript not found: {e}") from None
        if r.returncode != 0 or "NETMOSS_OK" not in r.stdout:
            raise RuntimeError(
                f"NetMoss2::getNetwork failed (exit {r.returncode})\n"
                f"{r.stderr[-4000:]}")
        with open(meta_fp) as f:
            run_meta = {k: v for k, v in json.load(f).items()}
        run_meta["input_sha256"] = fingerprint
        run_meta["wrapper_sha256"] = wrapper_fingerprint
        Path(meta_fp).write_text(json.dumps(run_meta, indent=2) + "\n")
    else:
        print(f"[{grp}] netmoss2 cache hit ({out_dir})")

    cons = pd.read_csv(pooled_fp, index_col=0)
    cons.index = cons.index.astype(str)
    cons.columns = cons.columns.astype(str)
    missing = set(fid) - set(cons.index)
    if missing:
        raise RuntimeError(f"pooled matrix missing {len(missing)} taxa, e.g. "
                           f"{sorted(missing)[:5]}")
    cons = cons.reindex(index=list(fid), columns=list(fid)).to_numpy()
    return cons, run_meta

def resolve_consensus(acc, units, fid, rscript, grp):
    """Consensus pooled matrix for sign reference (official NetMoss2 only).

    Non-finite NetMoss entries (saturated per-cohort r, no eps clip in the
    package) are patched with the in-house clipped value and counted.
    """
    cons_py = consensus_pool_py(acc)   # NaN-patch source (see its docstring)
    cons_r, meta = netmoss_pooled_matrix(units, fid, rscript, grp)
    engine = "netmoss2::getNetwork" \
        if meta.get("engine") == "netmoss2" else "netmoss::getNetwork"
    nf = ~np.isfinite(cons_r) & ~np.eye(len(cons_r), dtype=bool)
    n_sub = int(np.count_nonzero(nf))
    cons_r[nf] = cons_py[nf]          # in-house clipped value, mirrored side too
    off = ~np.eye(len(cons_r), dtype=bool)
    diff = cons_r[off] - cons_py[off]
    meta["py_vs_r_max_abs_diff"] = float(np.max(np.abs(diff))) if diff.size else 0.0
    meta["n_nonfinite_substituted_pairs"] = n_sub // 2
    meta.pop("n_substituted_nonfinite", None)
    print(f"[{grp}] netmoss vs in-house regression check: max|d|="
          f"{meta['py_vs_r_max_abs_diff']:.2e}, non-finite patched={n_sub // 2} pairs")
    return cons_r, engine, meta

def count_confirms(cohorts, F, cons_sign, r_thr):
    n_confirm = np.zeros((F, F), dtype=np.int32)
    for idx, r, w, n, u in cohorts:
        r_full = np.zeros((F, F))
        r_full[np.ix_(idx, idx)] = r
        strong = (np.abs(r_full) >= r_thr)
        confirm = strong & (np.sign(r_full) == cons_sign) & (cons_sign != 0)
        np.fill_diagonal(confirm, False)
        n_confirm += confirm
    return n_confirm

def build_tables(acc, cohorts, F, fid, r_thr, k_min, repro_thr, consensus_thr,
                 p0, cons_sign=None):
    if cons_sign is None:
        cons_sign = consensus_sign(acc)
    n_confirm = count_confirms(cohorts, F, cons_sign, r_thr)

    iu = np.triu_indices(F, k=1)
    cov_u = acc["cov"][iu]
    W1_u = acc["W1"][iu]
    VR_u = acc["VR"][iu]

    with np.errstate(invalid="ignore", divide="ignore"):
        pooled = np.where(W1_u > 0, VR_u / W1_u, np.nan)

    rmax = np.zeros((F, F))
    w_at = np.zeros((F, F))
    v_at = np.zeros((F, F))
    for idx, r, w, n, u in cohorts:
        ix = np.ix_(idx, idx)
        rabs = np.abs(r)
        blk = rmax[ix]
        newmax = rabs > blk
        rmax[ix] = np.where(newmax, rabs, blk)
        w_at[ix] = np.where(newmax, w, w_at[ix])
        v_at[ix] = np.where(newmax, r, v_at[ix])
    W_at_u = w_at[iu]
    V_at_u = v_at[iu]

    with np.errstate(invalid="ignore", divide="ignore"):
        pooled_wo = np.where(W1_u - W_at_u > 0,
                             (VR_u - W_at_u * V_at_u) / (W1_u - W_at_u), np.nan)

    nconf_u = n_confirm[iu]
    rep_u = np.divide(nconf_u, cov_u, out=np.full(cov_u.shape, np.nan), where=cov_u > 0)

    cand = cov_u >= k_min
    base = cand & (rep_u >= repro_thr) & (np.abs(pooled) >= consensus_thr)
    final = base

    binom_p = np.full(F * (F - 1) // 2, np.nan)
    binom_p[cand] = stats.binom.sf(nconf_u[cand] - 1, cov_u[cand], p0)
    binom_p_adj = np.full_like(binom_p, np.nan)
    if cand.any():
        binom_p_adj[cand] = bh_fdr(binom_p[cand])

    rows = np.flatnonzero(cand)
    ii, jj = iu[0][rows], iu[1][rows]
    df = pd.DataFrame({
        "taxon_i": fid[ii], "taxon_j": fid[jj],
        "pooled_obs": pooled[rows],
        "coverage": cov_u[rows].astype(int),
        "n_confirm": nconf_u[rows].astype(int),
        "repro": rep_u[rows],
        "binom_p": binom_p[rows],
        "binom_p_adj": binom_p_adj[rows],
        "top_weight_share": np.where(W1_u[rows] > 0,
                                     W_at_u[rows] / W1_u[rows], np.nan),
        "pooled_wo": pooled_wo[rows],
        "baseline": base[rows],
        "final": final[rows],
    })

    counts = dict(baseline=int(base.sum()), final=int(final.sum()))
    return df, counts, pooled, iu

def run_pooling():
    units_all = pd.read_csv(RESULT_DIR / "units.csv", dtype={"study": str})
    for group in GROUPS:
        units = discover_units(group, units_all)
        if not units:
            raise ValueError(f"No SparCC cohorts for {group}")
        features = collect_fid(units)
        accumulated, cohorts = build_stats(units, features)
        consensus, engine, metadata = resolve_consensus(accumulated, units, features, RSCRIPT, group)
        signs = np.sign(consensus)
        np.fill_diagonal(signs, 1)
        candidates, counts, _, _ = build_tables(accumulated, cohorts, len(features), features,
                                               0.20, 5, 0.60, 0.20, 0.025, signs)
        output = RESULT_DIR / "backbone_v2" / f"{group}_genus_uw"
        output.mkdir(parents=True, exist_ok=True)
        columns = ["taxon_i", "taxon_j", "pooled_obs", "coverage", "n_confirm", "repro",
                   "binom_p", "binom_p_adj", "top_weight_share", "pooled_wo"]
        candidates.loc[candidates.final, columns].to_csv(output / "backbone_edges.csv", index=False)
        parameters = dict(group=group, n_units=len(units), k_min=5, r_thr=0.20,
                          repro_thr=0.60, consensus_thr=0.20, p0=0.025,
                          pooling_engine=engine, pooling_detail=metadata, **counts)
        (output / "counts.json").write_text(json.dumps(parameters, indent=2) + "\n")
        print(f"[pool] {group}: {counts['final']} edges")
