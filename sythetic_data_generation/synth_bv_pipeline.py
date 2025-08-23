"""
Synthetic BV generator using ONLY user-input features (no file I/O).
Data will change depending on seed value.
Model:
  Y_isr = mu + B_i + W_is + A_isr
  B_i   ~ N(0, var_BP)         # between-subject
  W_is  ~ N(0, var_WP)         # within-subject (sample-level)
  A_isr ~ N(0, var_A)          # replicate/analytical

You provide:
  • mu, var_A, var_WP, var_BP  (from your original study)
  • a grid: balanced I×S×R   (or pass a custom grid DataFrame)
  • (optional) scale targets (mean/sd or min/max), rounding decimals, positivity clip

Outputs a pandas DataFrame with new synthetic "Result" values.
"""

from __future__ import annotations
import math
import numpy as np
import pandas as pd
from typing import Optional, Literal

# ----------------------------- grid helpers -----------------------------

def make_balanced_grid(I: int, S: int, R: int, subject_prefix: str = "P") -> pd.DataFrame:
    """Balanced Subject×Sample×Replicate grid."""
    rows = []
    for i in range(1, I + 1):
        subj = f"{subject_prefix}{i:02d}"
        for s in range(1, S + 1):
            for r in range(1, R + 1):
                rows.append((subj, s, r))
    return pd.DataFrame(rows, columns=["Subject", "Sample", "Replicate"])

def normalize_grid(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure grid has exactly ['Subject','Sample','Replicate'] with Subject=str, Sample=int, Replicate=int."""
    need = {"Subject","Sample","Replicate"}
    if not need.issubset(df.columns):
        raise KeyError("Custom grid must have columns: Subject, Sample, Replicate")
    out = df.loc[:, ["Subject","Sample","Replicate"]].copy()
    out["Subject"]   = out["Subject"].astype(str)
    out["Sample"]    = pd.to_numeric(out["Sample"], errors="raise").astype(int)
    out["Replicate"] = pd.to_numeric(out["Replicate"], errors="raise").astype(int)
    return out.reset_index(drop=True)

# ------------------------ variance tools (optional) ---------------------

def _rescale_to_variance(x: np.ndarray, target_var: float) -> np.ndarray:
    """Rescale centered vector x so sample variance ≈ target_var."""
    x = np.asarray(x, dtype=float)
    if target_var <= 0 or x.size < 2:
        return np.zeros_like(x)
    xc = x - x.mean()
    cur = float(xc.var(ddof=1)) if xc.size > 1 else 0.0
    if cur <= 0:
        return np.random.normal(0.0, math.sqrt(target_var), size=x.size)
    return xc * math.sqrt(max(target_var,0.0) / cur)

def estimate_components_unbalanced(df: pd.DataFrame):
    """
    Quick check of realized components from synthetic data (no files needed).
    Returns (grand, var_A, var_WP, var_BP).
    """
    df = df.rename(columns=str.title)
    subj_mean = df.groupby("Subject")["Result"].mean()
    samp_mean = df.groupby(["Subject","Sample"])["Result"].mean()
    grand     = df["Result"].mean()

    r_ij = df.groupby(["Subject","Sample"])["Replicate"].nunique()
    S_i  = df.groupby("Subject")["Sample"].nunique()

    r_bar = r_ij.mean()
    S_bar = S_i.mean()

    ss_bp = ((samp_mean.index.get_level_values(0).map(subj_mean) - grand)**2 * r_ij).groupby(level=0).sum().sum()
    ss_wp = ((samp_mean - samp_mean.index.get_level_values(0).map(subj_mean))**2 * r_ij).sum()
    ss_a  = ((df["Result"] - samp_mean.loc[list(zip(df["Subject"], df["Sample"]))].values)**2).sum()

    I    = max(subj_mean.size, 1)
    df_bp = max(I - 1, 1)
    df_wp = max((S_i - 1).sum(), 1)
    df_a  = max((r_ij - 1).sum(), 1)

    ms_bp = ss_bp / df_bp
    ms_wp = ss_wp / df_wp
    ms_a  = ss_a  / df_a

    var_A  = ms_a
    var_WP = max((ms_wp - ms_a) / max(r_bar, 1), 0.0)
    var_BP = max((ms_bp - ms_wp) / max(S_bar * r_bar, 1), 0.0)

    return float(grand), float(var_A), float(var_WP), float(var_BP)

# ----------------------------- core simulator -----------------------------

def simulate_bv_from_features(
    *,
    mu: float,
    var_A: float,
    var_WP: float,
    var_BP: float,
    grid: Optional[pd.DataFrame] = None,
    I: Optional[int] = None,
    S: Optional[int] = None,
    R: Optional[int] = None,
    seed: Optional[int] = 1234,
    enforce_exact_component_vars: bool = True,
    # scale/format knobs (all user-provided, no files)
    scale_mode: Literal["none","mean_sd","minmax"] = "none",
    ref_mean: Optional[float] = None,
    ref_sd:   Optional[float] = None,
    ref_min:  Optional[float] = None,
    ref_max:  Optional[float] = None,
    clip_nonnegative: bool = False,
    rounding_decimals: Optional[int] = None,
) -> pd.DataFrame:
    """
    Generate synthetic BV data from user-supplied features ONLY.

    - Provide either a custom 'grid' OR a balanced grid via I,S,R.
    - Optionally match scale to your known mean/sd or min/max, and rounding.
    """
    rng = np.random.default_rng(seed)

    # 1) Build/validate grid
    if grid is None:
        if None in (I,S,R):
            raise ValueError("Provide a custom 'grid' OR all of I, S, R for a balanced grid.")
        df = make_balanced_grid(I,S,R)
    else:
        df = normalize_grid(grid)

    # 2) Random effects
    subjects = df["Subject"].unique().tolist()
    ss_index = pd.MultiIndex.from_frame(df[["Subject","Sample"]]).drop_duplicates()

    B = pd.Series(rng.normal(0.0, math.sqrt(max(var_BP,0.0)), size=len(subjects)), index=subjects)
    W = pd.Series(rng.normal(0.0, math.sqrt(max(var_WP,0.0)), size=len(ss_index)), index=ss_index)
    A = pd.Series(rng.normal(0.0, math.sqrt(max(var_A,0.0)),  size=len(df)),        index=df.index)

    # 3) Optionally force realized component variances ≈ targets
    if enforce_exact_component_vars:
        if len(B)>1: B = pd.Series(_rescale_to_variance(B.values, var_BP), index=B.index)
        if len(W)>1: W = pd.Series(_rescale_to_variance(W.values, var_WP), index=W.index)
        if len(A)>1: A = pd.Series(_rescale_to_variance(A.values, var_A),  index=A.index)

    # 4) Compose Y = mu + B_i + W_is + A_isr
    y = (
        mu
        + df["Subject"].map(B).to_numpy()
        + W.loc[pd.MultiIndex.from_arrays([df["Subject"].values, df["Sample"].values])].to_numpy()
        + A.to_numpy()
    )
    out = df.copy()
    out["Result"] = y.astype(float)

    # 5) Optional scale matching (all inputs come from the user)
    if scale_mode == "mean_sd":
        if ref_mean is None or ref_sd is None:
            raise ValueError("scale_mode='mean_sd' requires ref_mean and ref_sd.")
        cur_m = float(out["Result"].mean())
        cur_s = float(out["Result"].std(ddof=1)) if len(out)>1 else 0.0
        if cur_s > 0:
            out["Result"] = (out["Result"] - cur_m) * (ref_sd/cur_s) + ref_mean
        else:
            out["Result"] = ref_mean
    elif scale_mode == "minmax":
        if ref_min is None or ref_max is None or ref_max <= ref_min:
            raise ValueError("scale_mode='minmax' requires valid ref_min < ref_max.")
        cur_min, cur_max = float(out["Result"].min()), float(out["Result"].max())
        if cur_max > cur_min:
            scale = (ref_max - ref_min) / (cur_max - cur_min)
            out["Result"] = (out["Result"] - cur_min) * scale + ref_min

    # 6) Optional clipping & rounding
    if clip_nonnegative:
        out["Result"] = out["Result"].clip(lower=0.0)
    if isinstance(rounding_decimals, int):
        out["Result"] = out["Result"].round(rounding_decimals)

    return out

# ----------------------------- example usage -----------------------------
if __name__ == "__main__":
    # === Paste YOUR features here (from your original study) ===
    # Variance components (same units as your data)
    MU     = 6.762916666666667    # grand mean
    VAR_A  = 0.15610416666666668  # replicate/analytical variance (σ²_A)
    VAR_WP = 0.6699826388888888   # within-subject variance (σ²_WP)
    VAR_BP = 0.8247793291962177   # between-subject variance (σ²_BP)

    # Grid describing your study layout (choose ONE of the two options):
    # Option A) Balanced grid — set these to match your study's counts:
    I, S, R = 30, 4, 2   # <<< REPLACE with your subject/sample/replicate counts

    # Option B) Custom grid — provide a DataFrame with these exact columns:
    # custom_grid = pd.DataFrame({
    #     "Subject":  ["P01","P01","P01","P02","P02","P02","P02"],
    #     "Sample":   [1,1,2,1,1,2,3],
    #     "Replicate":[1,2,1,1,2,1,1],
    # })

    # Optional scale controls (enter YOUR known scale, or leave off):
    #   - to match overall mean/sd:
    SCALE_MODE   = "mean_sd"      # "none", "mean_sd", or "minmax"
    REF_MEAN     = 6.762916666666667  # set to your dataset's mean if you want exact scale
    REF_SD       = 1.0             # set to your dataset's SD if known (example value)
    #   - OR to match min/max range:
    # SCALE_MODE = "minmax"; REF_MIN = 3.5; REF_MAX = 10.1

    # Generate synthetic data (no file reads)
    df_syn = simulate_bv_from_features(
        mu=MU, var_A=VAR_A, var_WP=VAR_WP, var_BP=VAR_BP,
        I=I, S=S, R=R,            # or grid=custom_grid
        seed=2025,
        enforce_exact_component_vars=True,
        scale_mode=SCALE_MODE,
        ref_mean=REF_MEAN, ref_sd=REF_SD,
        # ref_min=REF_MIN, ref_max=REF_MAX,
        clip_nonnegative=True,
        rounding_decimals=1,      # set to your typical decimals (e.g., 1)
    )

    print(df_syn.head())

    # (Optional) verify realized components on the synthetic data:
    gm, va, vwp, vbp = estimate_components_unbalanced(df_syn.rename(columns=str.title))
    print("\nRealized components on synthetic data:")
    print(f"  grand mean ≈ {gm:.6f}")
    print(f"  var_A      ≈ {va:.6f}")
    print(f"  var_WP     ≈ {vwp:.6f}")
    print(f"  var_BP     ≈ {vbp:.6f}")
    
df_syn
