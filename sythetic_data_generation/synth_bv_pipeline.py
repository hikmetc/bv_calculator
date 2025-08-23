# synth_bv_pipeline.py
# Reproduce the final synthetic dataset (see comments for steps).
# Requirements: pandas, numpy, openpyxl, xlsxwriter

import io
import os
import re
import numpy as np
import pandas as pd

# --------------------------- utilities ---------------------------

def normalize(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name).strip().lower())

def find_flexible(columns, candidates):
    cand_norm = [normalize(c) for c in candidates]
    for c in columns:                      # exact normalized match
        if normalize(c) in cand_norm:
            return c
    for c in columns:                      # substring match
        cn = normalize(c)
        if any(tok in cn for tok in cand_norm):
            return c
    return None

def detect_cols(df: pd.DataFrame):
    subj = find_flexible(df.columns, ["subject", "sunject", "participant", "patient", "id", "subjectid", "sunjectid"])
    samp = find_flexible(df.columns, ["sample", "time", "visit", "sampleid"])
    repl = find_flexible(df.columns, ["replicate", "rep", "replicateid", "repno"])
    res  = find_flexible(df.columns, ["result", "results", "value", "measurement"])
    if not all([subj, samp, repl, res]):
        missing = ["Subject","Sample","Replicate","Result"]
        have = {subj:"Subject", samp:"Sample", repl:"Replicate", res:"Result"}
        have_str = ", ".join([v for k,v in have.items() if k])
        raise KeyError(f"Could not detect all columns. Found: {have_str}. Needed: {', '.join(missing)}")
    df = df.rename(columns={subj:"Subject", samp:"Sample", repl:"Replicate", res:"Result"})
    return df[["Subject","Sample","Replicate","Result"]].copy()

def strip_ghost_index(df: pd.DataFrame) -> pd.DataFrame:
    if not df.empty and (df.columns[0] == "" or str(df.columns[0]).startswith("Unnamed")):
        return df.iloc[:,1:].copy()
    return df

def integer_like(series: pd.Series) -> bool:
    vals = pd.to_numeric(series.dropna(), errors="coerce")
    return vals.notna().all() and np.allclose(vals.values, np.rint(vals.values), atol=1e-9)

# ------------------- step 1: change per-group distributions -------------------

def change_results_per_group(df: pd.DataFrame, rng_seed: int = 12345) -> pd.Series:
    """
    Stratify by Subject×Sample×Replicate and replace only Result values:
      - Numeric: perturb mean & SD per group, sample Normal(new_mu, new_sigma),
                 clip to global range with ±10% margin.
      - Categorical: sample from a Dirichlet-perturbed category distribution.
    """
    np.random.seed(rng_seed)
    orig = df["Result"]
    mask = orig.isna()

    coerced = pd.to_numeric(orig, errors="coerce")
    numeric_share = coerced.notna().mean()

    global_nonnull = coerced[coerced.notna()]
    global_min = float(global_nonnull.min()) if not global_nonnull.empty else -np.inf
    global_max = float(global_nonnull.max()) if not global_nonnull.empty else  np.inf
    global_std = float(global_nonnull.std(ddof=0)) if not global_nonnull.empty else 1.0
    rng = global_max - global_min if np.isfinite(global_max) and np.isfinite(global_min) else None

    new_series = orig.copy()

    for _, idx in df.groupby(["Subject","Sample","Replicate"], dropna=False).indices.items():
        idx = pd.Index(idx)
        idx_nonnull = idx[~mask.loc[idx]]
        if len(idx_nonnull) == 0: 
            continue

        if numeric_share >= 0.8:
            g = coerced.loc[idx_nonnull]
            mu = float(g.mean()); sigma = float(g.std(ddof=0))
            base_std = sigma if (np.isfinite(sigma) and sigma > 0) else (global_std if (np.isfinite(global_std) and global_std > 0) else 1.0)
            shift = np.random.uniform(-0.75, 0.75) * base_std
            scale = np.random.uniform(0.6, 1.8)
            new_mu = mu + shift
            new_sigma = max(1e-9, base_std * scale)
            samples = np.random.normal(loc=new_mu, scale=new_sigma, size=len(idx_nonnull))
            if rng is not None and np.isfinite(rng) and rng > 0:
                lo = global_min - 0.1*rng; hi = global_max + 0.1*rng
                samples = np.clip(samples, lo, hi)
            new_series.loc[idx_nonnull] = samples
        else:
            g = orig.loc[idx_nonnull]
            counts = g.value_counts(dropna=True)
            alpha = counts.values.astype(float) + 0.5
            p = np.random.dirichlet(alpha)
            cats = counts.index.tolist()
            new_series.loc[idx_nonnull] = np.random.choice(cats, size=len(idx_nonnull), p=p)

    return new_series

# --------------------- step 2/3: rounding & variance shrink -------------------

def round_1dp(s: pd.Series) -> pd.Series:
    num = pd.to_numeric(s, errors="coerce").round(1)
    out = num.astype(float)
    # if original had non-numeric strings, preserve them
    nonnum_mask = pd.to_numeric(s, errors="coerce").isna() & s.notna()
    out.loc[nonnum_mask] = s.loc[nonnum_mask]
    return out

def shrink_within_subject_sample(df: pd.DataFrame, alpha: float) -> pd.Series:
    """
    Shrink deviations from the Subject×Sample mean by factor alpha (0<alpha<1):
       new = mean + alpha * (old - mean)
    """
    orig = df["Result"]
    num = pd.to_numeric(orig, errors="coerce")
    out = orig.copy()

    for (_, _), idx in df.groupby(["Subject","Sample"], dropna=False).indices.items():
        idx = pd.Index(idx)
        msk = num.loc[idx].notna()
        if not msk.any(): 
            continue
        sub_idx = idx[msk.values]
        gmean = num.loc[sub_idx].mean()
        out.loc[sub_idx] = gmean + alpha * (num.loc[sub_idx] - gmean)

    return round_1dp(out)

# --------------------------- main pipeline function ---------------------------

def generate_final_synthetic(in_path: str, out_path: str):
    # read CSV/XLSX
    _, ext = os.path.splitext(in_path.lower())
    if ext == ".xlsx":
        df = pd.read_excel(in_path)
    elif ext == ".csv":
        df = pd.read_csv(in_path)
    else:
        raise ValueError("Input must be .csv or .xlsx")

    df = strip_ghost_index(df)
    df = detect_cols(df)  # map to Subject, Sample, Replicate, Result

    # STEP 1: change per-group distributions (Subject×Sample×Replicate)
    df["Result"] = change_results_per_group(df, rng_seed=12345)

    # STEP 2: round to exactly 1 decimal
    df["Result"] = round_1dp(df["Result"])

    # STEP 3a: reduce replicate-based variance by ~30% (alpha=0.7)
    df["Result"] = shrink_within_subject_sample(df, alpha=0.7)  # keeps 1 dp

    # STEP 3b: additional 50% reduction (alpha=0.5) on the current values
    df["Result"] = shrink_within_subject_sample(df, alpha=0.5)  # keeps 1 dp

    # write Excel with column formatting = one decimal
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        df.to_excel(writer, sheet_name="Sheet1", index=False)
        ws = writer.sheets["Sheet1"]
        # apply number format to the Result column
        res_col_idx = df.columns.get_loc("Result")
        fmt = writer.book.add_format({"num_format": "0.0"})
        ws.set_column(res_col_idx, res_col_idx, None, fmt)

    return df

# ------------------------------ example usage --------------------------------
if __name__ == "__main__":
    # Change these paths as needed:
    input_file  = "template_2.xlsx"  # or your CSV/XLSX
    output_file = "results_replicate_variance_shrunk_1dp.xlsx"
    generate_final_synthetic(input_file, output_file)
    print(f"Saved: {output_file}")
