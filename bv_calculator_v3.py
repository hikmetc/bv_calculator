"""
Biological‑Variation Calculator (Streamlit)
==========================================
This single‑file Streamlit application helps a laboratory scientist calculate:

* **CV_A** – analytical imprecision
* **CV_I** – within‑subject biological variation (+95 % CI)
* **CV_G** – between‑subject biological variation
* **RCV** – two‑sided 95 % reference‑change value

It implements the *balanced* two‑level variance‑component algebra of **Røraas et al. (2012)**.

v5 highlights
-------------
* **Deep inline commentary** – every major line now has a `#` explanation so readers can
  follow both the Streamlit UI and the statistical maths.
* Template download in **CSV & XLSX** formats.
* Three data‑entry modes: *upload*, *paste*, *manual table*.

Run locally
-----------
```bash
pip install streamlit pandas numpy scipy openpyxl
streamlit run bv_streamlit.py
```
"""

# ——————————————————————————————————————————————————————————————
# 0.  Imports & general setup
# ——————————————————————————————————————————————————————————————
from __future__ import annotations   # future‑proof typing (PEP 563 style semantics)

import io                                # in‑memory byte streams for XLSX template
import textwrap                          # clean multi‑line strings for CSV template
from dataclasses import dataclass, field    # convenient result bundle
import numpy as np                       # scientific number‑crunching
import pandas as pd                      # tabular data handling
import streamlit as st                   # web UI framework
from scipy.stats import chi2             # chi‑square for exact CI of variance
from scipy.stats import shapiro, kstest, f              # normality + Cochran
from scipy.stats import bartlett, linregress  
import plotly.graph_objects as go   # ← NEW

# ——————————————————————————————————————————————————————————————
# 1.  Streamlit page config (title, layout)
# ——————————————————————————————————————————————————————————————
# Page configuration
st.set_page_config(
    page_title="Biological Variation Calculator",  # shown in browser tab
    layout="wide",
    initial_sidebar_state="expanded"                           # slim white gutter L/R
)

# ——————————————————————————————————————————————————————————————
# A.  Simple design system (tailwind-ish utility CSS)
# ——————————————————————————————————————————————————————————————
# Minimalistic CSS with enhanced heading
# Skeuomorphic CSS with white page background
st.markdown(
    """
    <style>
      /* Page background: plain white */
      body {
        background: #ffffff;
      }
      .block-container {
        padding: 2rem; 
        background: #ffffff;
        border: 1px solid #b0c4de;
        border-radius: 1rem;
        box-shadow: inset 0 1px 2px rgba(255,255,255,0.8), 0 4px 6px rgba(0,0,0,0.05);
        font-family: 'Georgia', serif;
        max-width: 900px;
        margin: 2rem auto;
      }
      /* Embossed heading in blue tones */
      h1 {
        text-align: center;
        font-size: 3rem;
        color: #1e3a5f;
        text-shadow: 1px 1px 2px rgba(255,255,255,0.7), -1px -1px 2px rgba(0,0,0,0.1);
        margin-bottom: 2rem;
      }
      /* Skeuomorphic metric cards with blue gradients */
      .metric-card {
        background: linear-gradient(145deg, #e8f4fd, #d0eaff);
        border-radius: 1rem;
        padding: 1.25rem;
        box-shadow: 4px 4px 8px rgba(0,0,0,0.1), inset 2px 2px 4px rgba(255,255,255,0.9);
        text-align: center;
        border: 1px solid #a3cce3;
      }
      .metric-title {
        font-size: 0.75rem;
        color: #28527a;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        margin-bottom: 0.5rem;
        text-shadow: 0 1px 1px rgba(255,255,255,0.9);
      }
      .metric-value {
        font-size: 2rem;
        color: #1b263b;
        font-weight: bold;
        text-shadow: 0 1px 1px rgba(255,255,255,0.8);
      }
      .metric-delta {
        font-size: 0.85rem;
        color: #3e5c76;
        margin-top: 0.5rem;
      }
      /* Button styles: soft blue sheen */
      .stButton>button {
        background: linear-gradient(145deg, #b8d6f6, #95c5f0);
        border-radius: 0.75rem;
        border: 1px solid #7fa6d8;
        padding: 0.75rem 1.5rem;
        box-shadow: 3px 3px 6px rgba(0,0,0,0.1), inset 2px 2px 4px rgba(255,255,255,0.8);
        font-family: 'Georgia', serif;
        color: #1e3a5f;
        font-weight: bold;
      }
      .stButton>button:hover {
        background: linear-gradient(145deg, #95c5f0, #b8d6f6);
      }
      /* DataFrame container */
      .stDataFrame {
        border: 1px solid #b0c4de;
        border-radius: 0.75rem;
        box-shadow: inset 0 1px 2px rgba(255,255,255,0.8);
      }
      .stDataFrame thead tr th {
        background: #d0e7f5;
        color: #1e3a5f;
        text-shadow: 0 1px 1px rgba(255,255,255,0.9);
      }
    </style>
    """,
    unsafe_allow_html=True
)

# ─── custom exception that stores the log ──────────────────────────────
class PreprocessError(ValueError):
    """Raised when the cleaner rejects the entire dataset.

    Attributes
    ----------
    log : list[str]
        The chronological quality-improvement log generated up to the
        failure point.  This lets the UI display what happened.
    """
    def __init__(self, msg: str, log: list[str]):
        super().__init__(msg)
        self.log = log

# ——————————————————————————————————————————————————————————————
# 2.  Dataclass to hold all computed results
# ——————————————————————————————————————————————————————————————
@dataclass
class BVResult:
    """Container for variance components and derived metrics."""

    # experiment counts
    I: int      # subjects
    S: int      # samples per subject
    R: int      # replicates per sample

    # mean of all raw results – needed for CV calculation
    grand_mean: float

    # variance components (sigma‑squared)
    var_A: float   # analytical
    var_WP: float  # within‑person
    var_BP: float  # between‑person

    # coefficients of variation (%)
    cv_A: float
    cv_I: float
    cv_G: float

    # confidence intervals & RCV
    ci_cv_I: tuple[float, float]  # (lower, upper) on CV_I
    ci_cv_G: tuple[float, float]  # (lower, upper) on CV_G   ← new
    rcv_95: float                 # % change considered significant (95 %)

    # NEW – verbose audit trail of what happened during cleaning
    preprocess_log: list[str] = field(default_factory=list)


# ——————————————————————————————————————————————————————————————
# 2-bis.  Column-mapping helpers
# ——————————————————————————————————————————————————————————————
_DEFAULT_GUESSES = {
    "Subject":   ("subject", "patient", "participant", "id", "birey"),
    "Sample":    ("sample", "time", "visit", "örnek"),
    "Replicate": ("rep", "replicate", "duplicate"),
    "Result":    ("result", "value", "measurement", "sonuç", "değer"),
}

def _default_idx(role: str, cols: list[str]) -> int:
    """Return index in *cols* whose cleaned name matches expected role."""
    keys = _DEFAULT_GUESSES[role]
    for i, c in enumerate(cols):
        if str(c).strip().lower().replace(" ", "") in keys:
            return i
    return 0


# ——————————————————————————————————————————————————————————————
# 3-bis.  Braga–Panteghini flow-chart utilities
#        (outlier tests + distribution checks)  :contentReference[oaicite:0]{index=0}
# ——————————————————————————————————————————————————————————————
def _cochrans_test(variances: np.ndarray,
                   alpha: float = 0.05,
                   df: int = 1) -> bool:
    """
    True  → largest variance is a Cochran outlier and must be rejected
    False → no Cochran outlier detected.
    """
    k = variances.size
    if k < 2:
        return False
    G = variances.max() / variances.sum()
    # critical limit from F-distribution (approximate Cochran)  :contentReference[oaicite:1]{index=1}
    q        = 1 - alpha / k
    f_crit   = f.ppf(q, 1, (k - 1) * df)
    G_crit   = f_crit / (f_crit + (k - 1))
    return G > G_crit


def _reed_outlier(values: np.ndarray) -> int | None:
    """
    Implements Reed’s criterion:
    returns index in *values* of the outlier, or None if none found.
    """
    if values.size < 3:
        return None
    s_idx = np.argsort(values)
    s_val = values[s_idx]
    rng   = s_val[-1] - s_val[0]
    if rng == 0:
        return None
    # high-end
    if (s_val[-1] - s_val[-2]) > rng / 3:
        return int(s_idx[-1])
    # low-end
    if (s_val[1] - s_val[0]) > rng / 3:
        return int(s_idx[0])
    return None


# ──────────────────────────────────────────────────────────────────────────────
# 3-bis.  Data cleaning & quality-improvement pipeline
#         (Braga-Panteghini flow-chart, QI 7-10, Reed, Cochran, etc.)
# ──────────────────────────────────────────────────────────────────────────────
def _preprocess_bv_dataframe(
        df_in: pd.DataFrame,
        alpha: float = 0.05,
        normal_p: float = 0.05
) -> tuple[pd.DataFrame, list[str]]:
    """
    Apply all pre-analysis quality-improvement (QI) checks and statistical
    filters required by Braga-Panteghini.  Returns a cleaned DataFrame **plus**
    a detailed chronological log (list of strings).
    """
    log: list[str] = []
    df = df_in.copy()

    # ════════════════════════════════════════════════════════════════════════
    # 0.  Duplicate-pair consistency (QI 8a)
    # ════════════════════════════════════════════════════════════════════════
    if df["Replicate"].nunique() == 2:                  # only meaningful for duplicates
        piv   = (df.pivot_table(index=["Subject", "Sample"],
                                columns="Replicate",
                                values="Result")
                   .dropna())                           # drop incomplete pairs
        delta = (piv[1] - piv[2]).abs()
        thr   = delta.mean() + 3 * delta.std(ddof=0)
        out_pairs = delta[delta > thr].index            # MultiIndex (Subject, Sample)

        for subj, samp in out_pairs:
            r_values = (piv.loc[(subj, samp)]
                           .sort_index()                # replicate 1 first
                           .to_dict())
            log.append(
                f"QI 8a – duplicate-pair outlier removed → "
                f"Subject {subj}, Sample {samp}: "
                + ", ".join(f"R{int(r)} = {v:.2f}" for r, v in r_values.items())
                + f"  (|Δ| = {abs(r_values[1]-r_values[2]):.2f} > {thr:.2f})"
            )

        # actually drop every offending pair
        if len(out_pairs):
            mask = df.set_index(["Subject", "Sample"]).index.isin(out_pairs)
            df   = df[~mask]
            log.append(f"Total duplicate-pair exclusions (QI 8a): {len(out_pairs)}")

    # ════════════════════════════════════════════════════════════════════════
    # 0 b.  Steady-state trend test (QI 7)
    # ════════════════════════════════════════════════════════════════════════
    drift_subj = []
    for subj, g in df.groupby("Subject"):
        if g["Sample"].nunique() > 2:                   # need ≥3 time-points
            slp, _, _, pval, _ = linregress(g["Sample"], g["Result"])
            if pval < 0.05:
                drift_subj.append((subj, slp, pval))

    for subj, slp, pval in drift_subj:
        log.append(
            f"QI 7 – temporal drift: Subject {subj} has significant trend "
            f"(slope = {slp:+.3g}, p = {pval:.3g}) – excluded."
        )
    if drift_subj:
        df = df[~df["Subject"].isin([s for s, *_ in drift_subj])]
        log.append(f"Total subjects excluded for drift (QI 7): {len(drift_subj)}")

    # ════════════════════════════════════════════════════════════════════════
    # 1.  Iterative outlier removal (Cochran then Reed)
    # ════════════════════════════════════════════════════════════════════════
    while True:
        # ── 1-a. Cochran variance test (within-subject)
        subj_var = df.groupby("Subject")["Result"].var(ddof=1).dropna()
        if subj_var.size >= 2 and _cochrans_test(subj_var.values, alpha):
            culprit = subj_var.idxmax()
            vals    = df.loc[df["Subject"] == culprit, "Result"]
            log.append(
                f"Cochran variance outlier → removed Subject {culprit} "
                f"(σ²_wp = {subj_var.max():.3f}; values: "
                + ", ".join(f"{v:.2f}" for v in vals) + ")"
            )
            df = df[~df["Subject"].eq(culprit)]
            continue

        # ── 1-b. Reed mean test (between-subject)
        subj_mean = df.groupby("Subject")["Result"].mean()
        ridx = _reed_outlier(subj_mean.values)
        if ridx is not None:
            culprit = subj_mean.index[ridx]
            vals    = df.loc[df["Subject"] == culprit, "Result"]
            log.append(
                f"Reed mean outlier → removed Subject {culprit} "
                f"(mean = {subj_mean.iloc[ridx]:.2f}; values: "
                + ", ".join(f"{v:.2f}" for v in vals) + ")"
            )
            df = df[~df["Subject"].eq(culprit)]
            continue
        break  # nothing else to drop

    # ════════════════════════════════════════════════════════════════════════
    # 2.  Normality assessment  (subject-level + subject-means)
    # ════════════════════════════════════════════════════════════════════════
    def _fraction_subjects_normal(d: pd.DataFrame) -> float:
        n = d["Subject"].nunique()
        if n < 3:            # too few subjects – force log-transform attempt
            return 0.0
        good = sum(
            shapiro(g["Result"])[1] > normal_p
            for _, g in d.groupby("Subject")
            if g["Result"].size >= 3
        )
        return good / n

    transformed = False
    if _fraction_subjects_normal(df) <= 0.5:
        df["Result"] = np.log(df["Result"])
        transformed  = True
        log.append("Natural-log transform applied – fewer than 50 % of subjects "
                   "passed Shapiro-Wilk on raw data.")

    # ── subject-mean normality
    means = df.groupby("Subject")["Result"].mean()
    if means.size >= 3:
        p_sw = shapiro(means)[1]
        if p_sw <= normal_p:
            z      = (means - means.mean()) / means.std(ddof=0)
            p_ks   = kstest(z, "norm")[1]
            if p_ks <= normal_p and not transformed:
                df["Result"] = np.log(df["Result"])
                transformed  = True
                log.append("Natural-log transform applied after subject-mean "
                           "normality failure.")
                means = df.groupby("Subject")["Result"].mean()
                if means.size >= 3 and shapiro(means)[1] <= normal_p:
                    raise ValueError("Normality could not be achieved – stopping.")
            elif p_ks <= normal_p:
                raise ValueError("Normality could not be achieved even after "
                                 "log-transform.")
    else:
        log.append("Subject-mean normality step skipped – fewer than three "
                   "subjects remained.")

    # ════════════════════════════════════════════════════════════════════════
    # 3.  Variance-homogeneity across subjects (Bartlett, QI 10)
    # ════════════════════════════════════════════════════════════════════════
    if df["Subject"].nunique() >= 2:
        bart_p = bartlett(*[g["Result"].values for _, g in df.groupby("Subject")])[1]
        msg    = ("heterogeneous" if bart_p < 0.05 else "homogeneous")
        log.append(f"Bartlett test (QI 10): p = {bart_p:.3g} → variances {msg}.")
    else:
        log.append("Bartlett test skipped – fewer than two subjects remain.")

    # ════════════════════════════════════════════════════════════════════════
    # 4.  Force a perfectly balanced design for the ANOVA
    # ════════════════════════════════════════════════════════════════════════
    # ════════════════════════════════════════════════════════════════════════
    # 4.  Force a perfectly balanced design (equal S & R)
    #     – required for the closed-form two-level ANOVA that follows
    # ════════════════════════════════════════════════════════════════════════

    # 4-a ▸ SUBJECT-LEVEL balance  ───────────────────────────────────────────
    samp_cnt  = df.groupby("Subject")["Sample"].nunique()          # how many samples per subject
    target_S  = samp_cnt.mode().iat[0] if not samp_cnt.empty else 0

    off_subj  = samp_cnt[samp_cnt != target_S]                     # subjects that deviate
    if not off_subj.empty:
        for subj, nS in off_subj.items():
            log.append(
                f"Balance check – Subject **{subj}** contributes "
                f"{nS}/{target_S} required samples; subject **excluded** to keep "
                "a fully crossed design."
            )
        log.append(
            f"{len(off_subj)} subject(s) dropped because they lacked the modal "
            f"sample count *S = {target_S}*."
        )

    # retain only subjects with the correct sample count
    df = df[df["Subject"].isin(samp_cnt[samp_cnt == target_S].index)]

    # 4-b ▸ SAMPLE-LEVEL balance  (replicate count)  ─────────────────────────
    rep_cnt  = (df.groupby(["Subject", "Sample"])["Replicate"]
                .nunique()
                .reset_index(name="n"))
    target_R = rep_cnt["n"].mode().iat[0] if not rep_cnt.empty else 0

    bad_pairs = rep_cnt[rep_cnt["n"] != target_R]
    if not bad_pairs.empty:
        for _, row in bad_pairs.iterrows():
            log.append(
                f"Balance check – Subject **{row.Subject}**, Sample **{row.Sample}** "
                f"has {row.n}/{target_R} replicate measurements; sample **removed**."
            )
        log.append(
            f"{len(bad_pairs)} sample(s) discarded to enforce a uniform replicate "
            f"count *R = {target_R}* across all remaining subjects."
        )

    if not bad_pairs.empty:
        bad_idx = df.set_index(["Subject", "Sample"]).index.isin(
            bad_pairs.set_index(["Subject", "Sample"]).index)
        df = df[~bad_idx]


    # 4-c.  final sanity check
    I_fin = df["Subject"].nunique()
    if df.empty:
        raise PreprocessError(
            "All data were excluded by quality checks – nothing left to analyse.",
            log,
        )
    S_fin = df.groupby("Subject")["Sample"].nunique().iloc[0]
    R_fin = (df.groupby(["Subject", "Sample"])["Replicate"].nunique().iloc[0])

    log.append(f"✅ Balanced data set ready: {I_fin} subjects × "
               f"{S_fin} samples × {R_fin} replicates retained for ANOVA.")

    # ════════════════════════════════════════════════════════════════════════
    # 5.  Return cleaned data (back-transformed if needed) + log
    # ════════════════════════════════════════════════════════════════════════
    if transformed:
        df["Result"] = np.exp(df["Result"])

    return df, log





# ——————————————————————————————————————————————————————————————
# 4.  Core calculation routine (closed‑form Røraas algebra)
# ——————————————————————————————————————————————————————————————

def calculate_bv(df: pd.DataFrame, alpha: float = 0.05) -> BVResult:
    """Compute balanced‑design variance components & CVs.

    Parameters
    ----------
    df : DataFrame
        Must contain at least the four columns: Subject, Sample, Replicate, Result.
        Extra columns are ignored. Capitalisation is normalised.
    alpha : float
        Significance level for the CI on CV_I (0.05 → 95 % CI).
    """

    
    # 3.1 normalise column names to Title‑Case so user can type "subject" etc.
    df = df.rename(str.title, axis=1)
    required = {"Subject", "Sample", "Replicate", "Result"}
    if not required.issubset(df.columns):
        # fail fast if critical columns missing
        raise KeyError(f"Missing required columns: {required - set(df.columns)}")
    
    # 3.1  apply Braga–Panteghini outlier/normality pipeline  
    try:
        df, pp_log = _preprocess_bv_dataframe(df, alpha=alpha)
    except PreprocessError as e:
        # Propagate the cleaner’s detailed log upward
        raise PreprocessError(str(e), e.log) from None

    # 3.2 derive counts (I, S, R) and assert balance
    I = df["Subject"].nunique()                                          # how many people
    S = df.groupby("Subject")["Sample"].nunique().iloc[0]               # samples per person
    R = df.groupby(["Subject", "Sample"])["Replicate"].nunique().iloc[0]  # replicates per sample

    # confirm every subject has S samples and every sample has R replicates
    if not (df.groupby("Subject")["Sample"].nunique() == S).all():
        raise ValueError("Unbalanced design: unequal samples per subject.")
    if not (
        df.groupby(["Subject", "Sample"])["Replicate"].nunique() == R
    ).all():
        raise ValueError("Unbalanced design: unequal replicates per sample.")

    # 3.3 overall mean – denominator for CV (%).
    grand = df["Result"].mean()

    # 3.4 means for each hierarchical level
    subj_mean = df.groupby("Subject")["Result"].mean()                   # \bar Y_{i..}
    samp_mean = df.groupby(["Subject", "Sample"])["Result"].mean()      # \bar Y_{is.}

    # 3.5 sums‑of‑squares according to EMS table (Røraas Eq 2)
    ss_bp = R * S * ((subj_mean - grand) ** 2).sum()                      # between‑person
    ss_wp = R * (
        (samp_mean - subj_mean.loc[samp_mean.index.get_level_values(0)].values) ** 2
    ).sum()  # within‑person (samples)
    ss_a = (
        (df["Result"] - samp_mean.loc[list(zip(df["Subject"], df["Sample"]))].values) ** 2
    ).sum()  # replicate analytical error

    # 3.6 convert SS → MS by dividing by appropriate degrees of freedom
    ms_bp = ss_bp / (I - 1)               # df = I‑1
    ms_wp = ss_wp / (I * (S - 1))         # df = I*(S‑1)
    ms_a = ss_a / (I * S * (R - 1))       # df = I*S*(R‑1)

    # 3.7 back‑solve variance components (closed‑form)
    var_A = ms_a                                          # σ²_A
    var_WP = max((ms_wp - ms_a) / R, 0.0)                 # σ²_WP (truncate <0 to 0)
    var_BP = max((ms_bp - ms_wp) / (S * R), 0.0)          # σ²_BP

    # 3.8 convert to CVs (%). grand mean in denominator.
    cv_A = np.sqrt(var_A) / grand * 100
    cv_I = np.sqrt(var_WP) / grand * 100
    cv_G = np.sqrt(var_BP) / grand * 100

    # 3.9 exact 95 % CI for CV_I using Burdick & Graybill’s method
    df_wp = I * (S - 1)

    # two‐tailed χ² quantiles
    chi2_lower_q = chi2.ppf(alpha/2, df_wp)       # e.g. 0.025 quantile
    chi2_upper_q = chi2.ppf(1 - alpha/2, df_wp)   # e.g. 0.975 quantile

    # limits on MS_WP
    ms_wp_lower = (df_wp * ms_wp) / chi2_upper_q
    ms_wp_upper = (df_wp * ms_wp) / chi2_lower_q

    # convert MS limits to σ²_WP limits: (MS_limit − σ²_A) ÷ R
    ci_var_low = max((ms_wp_lower - var_A) / R, 0.0)
    ci_var_up  = max((ms_wp_upper - var_A) / R, 0.0)

    # finally back-transform to %CV
    ci_cv_I = (
        np.sqrt(ci_var_low) / grand * 100,
        np.sqrt(ci_var_up)  / grand * 100,
    )

    # — 3.10 exact 95 % CI for CV_G (between-subject) ——
    df_bp = I - 1
    chi2_low_bp = chi2.ppf(alpha/2, df_bp)
    chi2_up_bp  = chi2.ppf(1 - alpha/2, df_bp)

    # MS_BP contains σ²_A + Rσ²_WP + SRσ²_BP; subtract lower levels first
    adj_ms_bp = ms_bp - ms_wp

    ms_bp_lower = (df_bp * adj_ms_bp) / chi2_up_bp
    ms_bp_upper = (df_bp * adj_ms_bp) / chi2_low_bp

    # convert MS limits to σ²_BP limits: MS_limit / (S*R)
    ci_var_bp_low = max(ms_bp_lower / (S * R), 0.0)
    ci_var_bp_up  = max(ms_bp_upper / (S * R), 0.0)

    ci_cv_G = (
        np.sqrt(ci_var_bp_low) / grand * 100,
        np.sqrt(ci_var_bp_up)  / grand * 100,
    )

    # 3.11 reference‑change value (two‑sided, 95 %)
    rcv = 1.96 * np.sqrt(2) * np.sqrt(cv_A**2 + cv_I**2)

    # bundle everything and return
    return BVResult(
        I, S, R, grand,
        var_A, var_WP, var_BP,
        cv_A, cv_I, cv_G,
        ci_cv_I, ci_cv_G,    # now passing both intervals
        rcv,
        preprocess_log=pp_log,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Beautiful per-subject mean ± range (min–max) plot
# ─────────────────────────────────────────────────────────────────────────────
def plot_subject_ranges(clean_df: pd.DataFrame) -> go.Figure:
    """
    Build a compact, minimalistic plot:

    • X-axis  → Subject IDs (only those that survived cleaning)  
    • Y-axis  → Result scale (mean •, min–max whisker)

    Returns
    -------
    plotly.graph_objects.Figure
    """
    # ── summarise data ──────────────────────────────────────────────────────
    agg = (clean_df.groupby("Subject")["Result"]
                  .agg(mean="mean", min="min", max="max")
                  .sort_index()
                  .reset_index())

    agg["err_low"]  = agg["mean"] - agg["min"]
    agg["err_high"] = agg["max"]  - agg["mean"]

    # create user-friendly x-labels (e.g. “Subject 7”)
    agg["x_label"] = agg["Subject"].astype(str).apply(lambda s: f"Subject {s}")

    # ── figure ──────────────────────────────────────────────────────────────
    fig = go.Figure(
        go.Scatter(
            x=agg["x_label"],
            y=agg["mean"],
            mode="markers",
            marker=dict(size=6, color="#2a3f5f"),
            error_y=dict(
                type="data",
                symmetric=False,
                array=agg["err_high"],
                arrayminus=agg["err_low"],
                color="#9ea4ae",
                thickness=1,
                width=0            # no caps → thin minimal whisker
            ),
            hovertemplate=(
                "%{x}<br>"
                "Mean %{y:.2f}<br>"
                "Range %{customdata[0]:.2f} – %{customdata[1]:.2f}"
                "<extra></extra>"
            ),
            customdata=agg[["min", "max"]],
        )
    )

    # ── styling ─────────────────────────────────────────────────────────────
    fig.update_layout(
        template="simple_white",
        height=300 if len(agg) <= 8 else 360,      # stays small
        margin=dict(l=50, r=10, t=30, b=60),
        xaxis=dict(
            title="",              # no title; labels self-explanatory
            tickangle=-45,         # slanted for readability
            tickfont=dict(size=11),
            showgrid=False
        ),
        yaxis=dict(
            title="Result",
            zeroline=False,
            showgrid=False,
            ticks="outside",
            ticklen=4,
            tickcolor="#bfc1c4"
        ),
        font=dict(size=11),
        title=dict(
            text="Per-subject mean ± range",
            x=0.5, xanchor="center",
            font=dict(size=14, color="#2a3f5f")
        ),
        hoverlabel=dict(bgcolor="white", font_size=11),
        showlegend=False
    )

    return fig

# ——————————————————————————————————————————————————————————————
# 5.  Prepare template datasets for users to download (CSV & XLSX)
# ——————————————————————————————————————————————————————————————

def get_template_df(
    num_subjects: int   = 40,
    num_samples:  int   = 5,
    num_replicates: int = 2,
) -> pd.DataFrame:
    """
    Generate a synthetic, yet perfectly balanced, study-design template.

    • Every subject has *num_samples* time-points / samples  
    • Every sample is analysed in duplicate (2 replicates)  
    • Results follow a simple deterministic pattern so the file is
      reproducible and easy to eyeball:  
          base = 100 + subject_index  
          result = base + sample_index + (replicate_index-1)
    """
    rows = []
    for s_idx in range(1, num_subjects + 1):
        subj_id = f"P{s_idx:02d}"          # P01, P02, …, P40
        base    = 100 + s_idx              # gives between-subject spread
        for sample in range(1, num_samples + 1):
            for rep in range(1, num_replicates + 1):
                result = base + sample + (rep - 1)   # tiny within-sample shift
                rows.append((subj_id, sample, rep, result))

    return pd.DataFrame(
        rows, columns=["Subject", "Sample", "Replicate", "Result"]
    )

# 4.1 build in‑memory CSV & XLSX once; reuse for all download clicks
_template_df = get_template_df()
_template_csv = _template_df.to_csv(index=False)           # plain text string
_template_xlsx_io = io.BytesIO()
with pd.ExcelWriter(_template_xlsx_io, engine="openpyxl") as xl_writer:
    _template_df.to_excel(xl_writer, sheet_name="Template", index=False)
_template_xlsx_io.seek(0)  # reset file pointer for download

# ——————————————————————————————————————————————————————————————
# 6.  Streamlit user interface
# ——————————————————————————————————————————————————————————————

# 6.1  headline
st.title("Biological Variation Calculator")

# 6.2  sidebar – quick help & template downloads
with st.sidebar:
    st.header("Quick start")
    st.write(
        "1. Upload or paste your balanced study file (CSV/XLSX)\n"
        "2. Map the four required columns\n"
        "3. Click on **Calculate** button"
    )
    st.divider()
    st.subheader("Need an example file?")
    with open('./template/bv_data_template.xlsx', "rb") as template_file:
        template_byte = template_file.read()
    # download template excel file
    st.download_button(label="⬇️ Click to Download Template File",
                        data=template_byte,
                        file_name="template.xlsx",
                        mime='application/octet-stream')
    st.write("---")
    st.info('*Developed by Hikmet Can Çubukçu, MD, PhD, MSc, EuSpLM* <hikmetcancubukcu@gmail.com>')



# 6.3 three tabs – three data‑entry modes
upload_tab, entry_tab = st.tabs(["Upload", "Manual Entry"])

user_df: pd.DataFrame | None = None  # will hold whichever dataframe the user supplies

# — Tab 1: Upload —
with upload_tab:
    up_file = st.file_uploader("Upload CSV or XLSX", type=["csv", "xlsx"])
    if up_file is not None:
        try:
            # decide reader based on extension
            if up_file.name.lower().endswith("xlsx"):
                user_df = pd.read_excel(up_file)
            else:
                user_df = pd.read_csv(up_file)
            st.success("File loaded ✓ – please review below.")
        except Exception as e:
            st.error(f"Load error: {e}")


# — Tab 2: Manual Entry —
with entry_tab:
    st.write("Double‑click to edit cells. Use the ➕ menu on the right to add rows.")

    # store editable df in session_state so edits persist across reruns
    if "manual_df" not in st.session_state:
        st.session_state.manual_df = _template_df.head(2)  # tiny starter grid

    manual_df = st.data_editor(
        st.session_state.manual_df,
        num_rows="dynamic",         # allow row addition
        use_container_width=True,
        key="editor",
    )
    # save edits back to session
    st.session_state.manual_df = manual_df

    # only accept manual_df as user_df if at least 4 fully populated rows exist
    if manual_df.dropna().shape[0] >= 4:
        user_df = manual_df.copy()

# ——————————————————————————————————————————————————————————————
# 7. Preview & calculate button
# ——————————————————————————————————————————————————————————————
if user_df is not None:
    st.subheader("Data preview")
    st.dataframe(user_df, use_container_width=True)


    # — Column mapping (collapsed until user opens it) —
    with st.expander("⇢  Map / confirm columns", expanded=False):
        cols = list(user_df.columns)
        with st.form("map_form"):
            c1, c2 = st.columns(2)
            with c1:
                subj_sel = st.selectbox("Subject identifier",  cols,
                                        index=_default_idx("Subject", cols))
                samp_sel = st.selectbox("Sample / time-point", cols,
                                        index=_default_idx("Sample", cols))
            with c2:
                repl_sel = st.selectbox("Replicate",          cols,
                                        index=_default_idx("Replicate", cols))
                res_sel  = st.selectbox("Result / value",     cols,
                                        index=_default_idx("Result", cols))

            confirmed = st.form_submit_button("Save mapping")

        if confirmed:
            st.session_state.mapped_df = user_df.rename(
                columns={
                    subj_sel: "Subject",
                    samp_sel: "Sample",
                    repl_sel: "Replicate",
                    res_sel:  "Result",
                }
            )
            st.success("Mapping saved – you can now calculate.")



    if st.button("Calculate", type="primary"):
        try:
            df_for_calc = st.session_state.get("mapped_df")
            if df_for_calc is None:
                st.warning("Please open “Map / confirm columns” and save a mapping first.")
                st.stop()                       # abort this run cleanly

            try:
                res = calculate_bv(df_for_calc)


                # Key metrics — big bold labels, smaller numbers
                st.subheader("Key metrics")
                cols = st.columns(4, gap="large")
                metrics = [
                    ("CV<sub>A</sub>",  f"{res.cv_A:.2f} %",                None),
                    ("CV<sub>I</sub>",  f"{res.cv_I:.2f} %",  f"{res.ci_cv_I[0]:.2f}–{res.ci_cv_I[1]:.2f}% CI"),
                    ("CV<sub>G</sub>",  f"{res.cv_G:.2f} %",  f"{res.ci_cv_G[0]:.2f}–{res.ci_cv_G[1]:.2f}% CI"),
                    ("95 % RCV",        f"±{res.rcv_95:.2f} %",             None),
                ]

                for col, (label, value, delta) in zip(cols, metrics):
                    with col:
                        st.markdown(f"""
                        <div style="
                            background: #f0f0f3;
                            border-radius: 1rem;
                            padding: 0.8rem 1rem;
                            text-align: center;
                            /* dual shadow for neumorphism */
                            box-shadow:
                            6px 6px 12px rgba(0, 0, 0, 0.12),
                            -6px -6px 12px rgba(255, 255, 255, 0.85);
                        ">
                        <!-- Label -->
                        <div style="
                            font-size: 1.4rem;
                            font-weight: 800;
                            color: #1e3a5f;
                            margin-bottom: 0.25rem;
                        ">{label}</div>
                        <!-- Value -->
                        <div style="
                            font-size: 1rem;
                            font-weight: 600;
                            color: #1b263b;
                            line-height: 1.2;
                        ">{value}</div>
                        <!-- Delta/CI -->
                        {f"<div style='font-size:0.75rem; color:#5d6d7e; margin-top:0.3rem;'>{delta}</div>" if delta else ""}
                        </div>
                        """, unsafe_allow_html=True)


                # Summary table of CVs and 95% CIs
                st.subheader("Summary of variation metrics")
                var_tbl = pd.DataFrame({
                    "Component": ["Analytical", "Within-subject", "Between-subject"],
                    "Variance":  [res.var_A, res.var_WP, res.var_BP],
                    "CV %":      [res.cv_A, res.cv_I, res.cv_G],
                    "95 % CI":   ["–",
                                f"{res.ci_cv_I[0]:.2f}–{res.ci_cv_I[1]:.2f} %",
                                f"{res.ci_cv_G[0]:.2f}–{res.ci_cv_G[1]:.2f} %"],
                })
                st.dataframe(
                    var_tbl.style
                        .format({"Variance": "{:.3f}", "CV %": "{:.2f}"})
                        .set_table_styles([{
                            "selector": "th",
                            "props": [("background-color", "#d0e7f5"),
                                        ("color", "#1e3a5f"),
                                        ("font-size", "0.9rem")]
                        }]),
                    hide_index=True,
                    use_container_width=True
                )



                # — Per-subject mean ± range plot ————————————————————————————————
                clean_df, _ = _preprocess_bv_dataframe(df_for_calc)
                st.subheader("Per-subject distribution")
                st.plotly_chart(plot_subject_ranges(clean_df), use_container_width=True)





                # ─────────────────────────────────────────────────────────────────────────────
                #  📋  BIVAC CHECKLIST  (QI-6 → QI-13)
                #      Grades:  A (best) → D (unreliable)
                # ─────────────────────────────────────────────────────────────────────────────
                def _build_qi_checklist(log_lines: list[str], res: BVResult) -> pd.DataFrame:
                    """
                    Create a BIVAC quality table that covers QI-6 … QI-13.

                    • Every row starts at grade “A”.  
                    • We **only** downgrade when the *requirement itself* is not fulfilled –
                    finding / removing drifts or outliers does **NOT** lower the score
                    (BIVAC rewards the presence of those checks).  
                    • The “Details” column is populated with quantitative numbers and the
                    exact log lines that justify the given grade.
                    """

                    # ── canonical explanations (taken from the BIVAC PDF) ──────────────────
                    expl = {
                        "6": "Analytical imprecision documented (CVₐ from duplicates / IQC)",
                        "7": "Subjects in steady-state (individual time-trend examined)",
                        "8": "Comprehensive outlier handling (replicates, samples, subjects)",
                        "9": "Normal distribution achieved (raw or after transformation)",
                        "10": "Homogeneity of within-subject variances (Bartlett, p > 0.05)",
                        "11": "Appropriate 2-level (nested) ANOVA model applied",
                        "12": "Exact 95 % confidence limits reported",
                        "13": "Transparency on number of results retained vs. excluded",
                    }

                    # ── helper: make one blank† row  ────────────────────────────────────────
                    def _row(qi: str) -> dict[str, str]:
                        return dict(
                            QI=f"QI {qi}",
                            Grade="A",
                            Comment="Meets recommendation",
                            Explanation=expl[qi],
                            Details="",
                        )

                    rows = {q: _row(q) for q in expl}

                    # pre-populate quantitative fields (always valuable)
                    rows["6"]["Details"]  = f"CVₐ = {res.cv_A:.2f} %"
                    rows["12"]["Details"] = (
                        f"CVᵢ 95 % CI {res.ci_cv_I[0]:.2f}–{res.ci_cv_I[1]:.2f} %; "
                        f"CVg 95 % CI {res.ci_cv_G[0]:.2f}–{res.ci_cv_G[1]:.2f} %"
                    )
                    rows["13"]["Details"] = (
                        f"Subjects = {res.I} • "
                        f"Samples/subject = {res.S} • "
                        f"Replicates/sample = {res.R}  →  N = {res.I*res.S*res.R}"
                    )

                    # ── scan the verbose log to adjust grades / add notes ──────────────────
                    for raw in log_lines:
                        l = raw.lower()

                        # QI-7  (steady-state)
                        if "temporal drift" in l:
                            rows["7"]["Comment"] = "Drift detected; affected subject removed"
                            rows["7"]["Details"] += (" · " + raw)

                        # QI-8  (outliers)
                        if "duplicate-pair outlier" in l:
                            rows["8"]["Details"] += (" · " + raw)
                        if ("reed mean outlier" in l) or ("cochran variance outlier" in l):
                            rows["8"]["Details"] += (" · " + raw)

                        # QI-9  (normality)
                        if "log transform applied" in l:
                            rows["9"]["Comment"] = "Log-transform applied (allowed by BIVAC)"
                            rows["9"]["Details"] += (" · " + raw)
                        if "normality could not be achieved" in l:
                            rows["9"]["Grade"]   = "D"
                            rows["9"]["Comment"] = "Gaussian assumption ultimately failed"
                            rows["9"]["Details"] += (" · " + raw)

                        # QI-10  (variance homogeneity)
                        if "bartlett test" in l:
                            rows["10"]["Details"] += (" · " + raw)
                            if "heterogeneous" in l:
                                rows["10"]["Grade"]   = "C"
                                rows["10"]["Comment"] = "Variances heterogeneous (p < 0.05)"

                        # QI-13  (result transparency)
                        if ("subject excluded" in l) or ("sample removed" in l):
                            rows["13"]["Details"] += (" · " + raw)

                    # ── turn into tidy DataFrame – order by QI number, A→D categorical grade
                    df = (pd.DataFrame(rows.values())
                            .assign(Grade=lambda d: pd.Categorical(d["Grade"],
                                                                categories=list("ABCD"),
                                                                ordered=True))
                            .sort_values("QI")
                            .reset_index(drop=True))
                    return df




                # ── render the checklist & XLSX download ────────────────────────────────────
                st.subheader("BIVAC Checklist (QI 6 → QI 13)")
                bivac_df = _build_qi_checklist(res.preprocess_log, res)

                st.dataframe(
                    bivac_df.style.apply(
                        lambda s: ["background:#e8f9f0" if g in ("A", "B") else
                                "background:#fdecea"   for g in s],
                        axis=1, subset=["Grade"]
                    ).set_properties(**{"font-size": "0.85rem"}),
                    use_container_width=True,
                    hide_index=True,
                )

                # ⬇️  Excel export (tries xlsxwriter → falls back to openpyxl)
                with st.expander("⇢ Download checklist"):
                    import io, importlib
                    engine = "xlsxwriter" if importlib.util.find_spec("xlsxwriter") else "openpyxl"

                    xio = io.BytesIO()
                    with pd.ExcelWriter(xio, engine=engine) as xl:
                        bivac_df.to_excel(xl, index=False, sheet_name="BIVAC_QI_6-13")

                    st.download_button(
                        "Excel file",
                        data=xio.getvalue(),
                        file_name="BIVAC_QI_6-13.xlsx",
                        mime=(
                            "application/vnd.openxmlformats-officedocument."
                            "spreadsheetml.sheet"
                        ),
                    )





                # — Outlier / normality audit trail —
                with st.expander("Pre-processing log"):
                    if res.preprocess_log:
                        for line in res.preprocess_log:
                            st.write("•", line)
                    else:
                        st.write("No outliers detected; data met normality assumptions.")

            except PreprocessError as e:
                st.error(e.args[0])                      # friendly headline
                st.subheader("Quality-Improvement log (all steps)")
                for line in e.log:
                    st.write("•", line)                  # bullet list
            except Exception as e:
                st.error(f"Calculation failed: {e}")     # any other unexpected error

        except Exception as e:
            st.write(e)
            st.error(f"Calculation failed: {e}")
else:
    st.info("Input data above to enable calculation.")
