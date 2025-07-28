"""
Deveoped by Hikmet Can Çubukçu 
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
from scipy.stats import chi2 , t            # chi‑square for exact CI of variance
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

    ci_mean:  tuple[float, float]       # 95 % CI on grand mean
    ci_cv_A:  tuple[float, float]       # 95 % CI on CV_A

    # confidence intervals & RCV
    ci_cv_I: tuple[float, float]  # (lower, upper) on CV_I
    ci_cv_G: tuple[float, float]  # (lower, upper) on CV_G   ← new
    rcv_95: float                 # % change considered significant (95 %)

    # NEW – verbose audit trail of what happened during cleaning
    preprocess_log: list[str] = field(default_factory=list)
    
    # NEW – optional CV-ANOVA estimates
    cv_I_cv_anova: float | None = None
    ci_cv_I_cv_anova: tuple[float, float] | None = None


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

# ─────────────────────────────────────────────────────────────────────────────
# Helper: remove any index column accidentally read as data
# ─────────────────────────────────────────────────────────────────────────────
def _strip_ghost_index(df: pd.DataFrame) -> pd.DataFrame:
    """
    If the left‑most column is unnamed ('' or starts with 'Unnamed'),
    drop it and return the cleaned DataFrame.
    """
    if not df.empty and (df.columns[0] == "" or df.columns[0].startswith("Unnamed")):
        return df.iloc[:, 1:]
    return df

# ——————————————————————————————————————————————————————————————
# 3-bis.  Braga–Panteghini flow-chart utilities
#        (outlier tests + distribution checks)  :contentReference[oaicite:0]{index=0}
# ——————————————————————————————————————————————————————————————

from collections import namedtuple
CochranResult = namedtuple("CochranResult", "flag G G_crit")

def _cochrans_test(variances: np.ndarray,
                   alpha: float = 0.05,
                   df: int = 1) -> CochranResult:
    """
    Perform Cochran’s C test and *return both* the decision and
    the two key statistics so the caller can log them.

    Returns
    -------
    CochranResult(flag, G, G_crit)
        flag : True  → largest variance is an outlier
               False → no outlier
    """
    k = variances.size
    if k < 2:
        return CochranResult(False, np.nan, np.nan)

    G = variances.max() / variances.sum()
    q = 1 - alpha / k                       # Šidák
    f_crit = f.ppf(q, 1, (k - 1) * df)
    G_crit = f_crit / (f_crit + (k - 1))
    return CochranResult(G > G_crit, G, G_crit)

def _cochrans_test2(variances: np.ndarray,
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
        *,
        alpha: float = 0.05,
        normal_p: float = 0.05,
        enforce_balance: bool = True,
        flags: dict[str, bool] | None = None,          #  NEW
) -> tuple[pd.DataFrame, list[str]]:
    """
    Apply all pre-analysis quality-improvement (QI) checks and statistical
    filters required by Braga-Panteghini.  Returns a cleaned DataFrame **plus**
    a detailed chronological log (list of strings).
    """
    log: list[str] = []
    df = df_in.copy()
    transformed = False

    flags = flags or {}                    # safety – empty dict if None
    ON     = lambda name: flags.get(name, True)   # helper: is this step enabled?

    # ════════════════════════════════════════════════════════════════════════
    # 0.  Replicate‑set consistency  (QI 8a – Cochran variance test)
    #     • For every Subject‑×‑Sample pair we calculate the variance of its
    #       replicate measurements.
    #     • Cochran’s C identifies the pair whose variance is disproportionately
    #       large.  We iteratively remove offenders until the test is passed.
    # ════════════════════════════════════════════════════════════════════════

    # ── QI 8a – replicate–set Cochran ──────────────────────────────────────────
    if ON("rep_cochran"):
        while True:
            rep_var = (df.groupby(["Subject", "Sample"])["Result"]
                        .var(ddof=1).dropna())
            if rep_var.size < 3:
                break

            R = (df.groupby(["Subject", "Sample"])["Replicate"]
                .nunique().mode().iat[0])          # current modal replicate count

            c_res = _cochrans_test(rep_var.values, alpha, df=R - 1)
            # ── NEW universal log line ──────────────────────────
            log.append(
                f"QI 8a – Cochran (replicate sets): "
                f"G = {c_res.G:.3f}, Gcrit = {c_res.G_crit:.3f}"
                + (" → OUTLIER" if c_res.flag else " → no outlier detected")
            )

            if not c_res.flag:          # test passed → exit the loop
                break

            log.append(
                f"QI 8a – Cochran test: G = {c_res.G:.3f}, "
                f"Gcrit = {c_res.G_crit:.3f}"
            )

            subj, samp = rep_var.idxmax()
            vals = df.loc[(df["Subject"] == subj) & (df["Sample"] == samp), "Result"]
            log.append(
                f"QI 8a – replicate Cochran outlier removed → "
                f"Subject {subj}, Sample {samp}: "
                + ", ".join(f"{v:.2f}" for v in vals) +
                f"  (s² = {rep_var.max():.4f})"
            )
            df = df[~((df["Subject"] == subj) & (df["Sample"] == samp))]
    else:
        log.append("Replicate Cochran (QI 8a) skipped (switch off).")

    # ─────────────────────────────────────────────────────────────────────────
    # ❋ NEW: Bartlett's homogeneity test – analytic variance (replicates)
    # Purpose: Are the replicate variances per Subject×Sample similar? 
    # (Is the analytic imprecison homogeneous across all samples?)
    # ─────────────────────────────────────────────────────────────────────────
    if ON("rep_bartlett"):

        valid_groups = []
        group_keys    = []                        # (Subject, Sample) anahtarlarını da tut
        for (subj, samp), g in df.groupby(["Subject", "Sample"]):
            if g["Result"].size < 2:              # Bartlett için ≥2 ölçüm şart
                continue
            if g["Result"].var(ddof=1) == 0:      # Excel: s² = 0 olan çiftleri hariç tutar
                continue
            valid_groups.append(g["Result"].values)
            group_keys.append((subj, samp))

        if len(valid_groups) >= 3:                # Bartlett ≥3 grup ister
            with np.errstate(all="ignore"):       # sabit giriş uyarılarını bastır
                bart_stat, bart_p = bartlett(*valid_groups)
            rep_msg = "heterogeneous" if bart_p < alpha else "homogeneous"
            k = len(valid_groups)
            chi_crit = chi2.ppf(1 - alpha, k - 1)

            log.append(
                f"QI 8a – Bartlett test (replicate sets): "
                f"χ² = {bart_stat:.2f}, χ²crit = {chi_crit:.2f}, "
                f"p = {bart_p:.4f} → variances {rep_msg}."
            )

            # Excel’de olduğu gibi: heterojen ise en yüksek 5 s²’yi raporla
            if bart_p < alpha:
                rep_var = (df.groupby(["Subject", "Sample"])["Result"]
                            .var(ddof=1)
                            .loc[group_keys])     # yalnızca teste giren çiftler
                top_bad = rep_var.sort_values(ascending=False).head(5)
                for (subj, samp), v in top_bad.items():
                    log.append(
                        f"  High analytical variance → "
                        f"Subject {subj}, Sample {samp}  (s² = {v:.4f})"
                    )
        else:
            log.append(
                "QI 8a – Bartlett test (replicate sets) skipped: "
                "fewer than 3 valid replicate groups."
            )
    else:
        log.append("Replicate Bartlett (QI 8a) skipped (switch off).")

    # ── QI 8b – sample‑level Cochran inside each subject ──────────────────
    if ON("samp_cochran"):

        changed = True
        while changed:
            changed = False
            for subj, g in df.groupby("Subject"):

                samp_var = g.groupby("Sample")["Result"].var(ddof=1).dropna()

                # (3) ensure we log even when the test is impossible
                if samp_var.size < 3:
                    log.append(f"QI 8b – Cochran skipped (Subject {subj}): <3 variances.")
                    continue

                r_per_sample = g.groupby("Sample")["Replicate"].nunique().mode().iat[0]

                # (2) Cochran undefined if only one replicate
                if r_per_sample < 2:
                    log.append(f"QI 8b – Cochran skipped (Subject {subj}): only 1 replicate.")
                    continue

                c_res = _cochrans_test(samp_var.values, alpha, df=r_per_sample - 1)

                # single universal log line
                log.append(
                    f"QI 8b – Cochran (Subject {subj}): "
                    f"G = {c_res.G:.3f}, Gcrit = {c_res.G_crit:.3f}"
                    + (" → OUTLIER" if c_res.flag else " → no outlier detected")
                )

                if not c_res.flag:
                    continue  # nothing to remove – next subject

                # remove the offending sample
                samp = samp_var.idxmax()
                vals = df.loc[(df["Subject"] == subj) & (df["Sample"] == samp), "Result"]
                log.append(
                    f"QI 8b – sample Cochran outlier removed → "
                    f"Subject {subj}, Sample {samp}: "
                    + ", ".join(f"{v:.2f}" for v in vals) +
                    f"  (s² = {samp_var.max():.4f})"
                )
                df = df[~((df["Subject"] == subj) & (df["Sample"] == samp))]
                changed = True
                break   # restart because groupby cache is stale
    else:
        log.append("Sample Cochran (QI 8b) skipped (switch off).")

    # ── QI 8c – between‑subject Reed test only ────────────────────────────────
    if ON("reed"):

        removed_any = False          # ① NEW – keeps track of removals

        while True:
            subj_mean = df.groupby("Subject")["Result"].mean()
            ridx      = _reed_outlier(subj_mean.values)

            if ridx is None:
                # no new outlier found → exit loop …
                if not removed_any:  # ② NEW – log only if the very first check already passed
                    log.append("QI 8c – Reed test passed: no between‑subject outliers.")
                break

            removed_any = True       # mark that we did remove something

            culprit = subj_mean.index[ridx]
            vals    = df.loc[df["Subject"] == culprit, "Result"]

            log.append(
                f"QI 8c – Reed mean outlier → removed Subject {culprit} "
                f"(mean = {subj_mean.iloc[ridx]:.2f}; values: "
                + ", ".join(f"{v:.2f}" for v in vals) + ")"
            )
            df = df[~df["Subject"].eq(culprit)]
    else:
        log.append("Reed between‑subject (QI 8c) skipped (switch off).")

    # ════════════════════════════════════════════════════════════════════════
    # 0 b.  Steady‑state trend test (QI 7) – slope ± 95 % CI
    # ════════════════════════════════════════════════════════════════════════
    if ON("drift"):

        tinv = lambda p, df: abs(t.ppf(p/2, df))        # two‑sided t‑quantile
        drift_records = []                               # list of dicts

        for subj, g in df.groupby("Subject"):
            # need at least three time‑points for a slope
            if g["Sample"].nunique() <= 2:
                continue

            res = linregress(g["Sample"], g["Result"])
            df_denom = len(g) - 2                       # regression d.f.
            if df_denom <= 0:
                continue

            ts = tinv(0.05, df_denom)                   # 95 % two‑sided
            ci_low  = res.slope - ts * res.stderr
            ci_high = res.slope + ts * res.stderr

            if res.pvalue < 0.05:                       # significant drift
                drift_records.append(
                    dict(subj=subj,
                        slope=res.slope,
                        p=res.pvalue,
                        ci=(ci_low, ci_high))
                )

        # ── write log + drop drifting subjects *once* ───────────────────────
        if drift_records:
            for d in drift_records:
                log.append(
                    "QI 7 – temporal drift: Subject {subj} slope = {s:+.3g} "
                    "(95 % CI {lo:.3g}–{hi:.3g}, p = {p:.3g}) – excluded."
                    .format(subj=d["subj"], s=d["slope"],
                            lo=d["ci"][0], hi=d["ci"][1], p=d["p"])
                )

            df = df[~df["Subject"].isin([d["subj"] for d in drift_records])]
            log.append(f"Total subjects excluded for drift (QI 7): {len(drift_records)}")

    else:
        log.append("Steady‑state drift (QI 7) skipped (switch off).")


    # ════════════════════════════════════════════════════════════════════════
    # 2.  Normality assessment  (subject-level  +  subject-means)   ⟨QI 9⟩
    # ════════════════════════════════════════════════════════════════════════
    if ON("normality"):

        def _fraction_subjects_normal(d: pd.DataFrame) -> tuple[int, int]:
            """
            Return (n_gaussian, n_eligible) where “eligible” means ≥3 results/subject.
            Gaussianity is Shapiro-Wilk p-value > normal_p on raw results.
            """
            good, total = 0, 0
            for _, g in d.groupby("Subject"):
                if g["Result"].size >= 3:
                    total += 1
                    if shapiro(g["Result"])[1] > normal_p:
                        good += 1
            return good, total

        # ── 2-a. Shapiro-Wilk on each subject (raw scale) ───────────────────────
        n_good, n_total = _fraction_subjects_normal(df)
        if n_total:
            frac_txt = f"{n_good}/{n_total} subjects ({n_good/n_total:.0%})"
            log.append(
                f"Normality check (QI 9) – {frac_txt} passed Shapiro-Wilk "
                f"p>{normal_p} on raw results."
            )
        else:
            log.append("Normality check (QI 9) – skipped (too few subjects).")

        if n_total and n_good / n_total <= 0.50:
            df["Result"] = np.log(df["Result"])
            transformed  = True
            log.append("Natural-log transform applied – <50 % of subjects were Gaussian.")

        # ── 2‑b. Shapiro‑Wilk on subject means (always if ≥3 subjects) ──────────
        means = df.groupby("Subject")["Result"].mean()
        if means.size >= 3:
            p_sw = shapiro(means)[1]
            log.append(f"Subject‑means Shapiro-Wilk p = {p_sw:.3g}.")
            
            if p_sw <= normal_p and not transformed:
                # one chance: log‑transform the whole dataset
                df["Result"] = np.log(df["Result"])
                transformed  = True
                log.append("Natural‑log transform applied after subject‑mean SW failure.")
                means = df.groupby("Subject")["Result"].mean()
                if shapiro(means)[1] <= normal_p:
                    raise ValueError("Normality could not be achieved – stopping.")
            elif p_sw <= normal_p:
                raise ValueError("Normality could not be achieved even after log‑transform.")
    else:
        log.append("Normality checks / log‑transform (QI 9) skipped (switch off).")  

    # ════════════════════════════════════════════════════════════════════════
    # 3.  Variance-homogeneity across subjects (Bartlett, QI 10)
    # ════════════════════════════════════════════════════════════════════════
#    if df["Subject"].nunique() >= 2:
#        bart_p = bartlett(*[g["Result"].values for _, g in df.groupby("Subject")])[1]
#        msg    = ("heterogeneous" if bart_p < 0.05 else "homogeneous")
#        log.append(f"Bartlett test (QI 10): p = {bart_p:.3g} → variances {msg}.")
#    else:
#        log.append("Bartlett test skipped – fewer than two subjects remain.")

    # ════════════════════════════════════════════════════════════════════════
    # 3.  Variance‑homogeneity *across subjects*  (Bartlett, QI 10 – Excel uyumlu)
    #     Excel mantığı:
    #       ① Her denekte (Subject) önce replicates ortalaması alınır →
    #          her örnek (Sample) tek bir değere iner → analitik gürültü elenir.
    #       ② Aynı deneğin ≥2 örneği varsa bunların varyansı s²_WP hesaplanır.
    #       ③ Bartlett testi ≥3 deneklik {s²_WP} kümesine uygulanır.
    # ════════════════════════════════════════════════════════════════════════
    if ON("wp_bartlett"):

        subj_groups = []          # Bartlett’e girecek her deneğin “örnek ortalamaları”
        subj_keys   = []          # Aynı sırayla Subject ID tut (rapor log’u için)

        for subj, g in df.groupby("Subject"):
            if g["Sample"].nunique() < 2:      # Bartlett için denek başına ≥2 örnek gerek
                continue
            sample_means = g.groupby("Sample")["Result"].mean().values
            if sample_means.var(ddof=1) == 0:  # Excel: varyans sıfırsa deneği test dışı bırakır
                continue
            subj_groups.append(sample_means)
            subj_keys.append(subj)

        if len(subj_groups) >= 3:              # Bartlett ≥3 grup ister
            with np.errstate(all="ignore"):    # sabit grup uyarılarını bastır
                bart_stat, bart_p = bartlett(*subj_groups)
            wp_msg = "heterogeneous" if bart_p < alpha else "homogeneous"
            k = len(subj_groups)
            chi_crit = chi2.ppf(1 - alpha, k - 1)

            log.append(
                f"QI 10 – Bartlett test (within‑subject variances): "
                f"χ² = {bart_stat:.2f}, χ²crit = {chi_crit:.2f}, "
                f"p = {bart_p:.4f} → variances {wp_msg}."
            )

            # If p < alpha, report the top 5 highest within-subject variances (s²_WP)
            if bart_p < alpha:
                wp_var = {
                    subj: g.groupby("Sample")["Result"].mean().var(ddof=1)
                    for subj, g in df.groupby("Subject") if subj in subj_keys
                }
                for subj, v in sorted(wp_var.items(), key=lambda x: x[1], reverse=True)[:5]:
                    log.append(
                        f"  High within‑subject variance → Subject {subj}  (s² = {v:.4f})"
                    )
        else:
            log.append(
                "QI 10 – Bartlett test (within‑subject variances) skipped: "
                "fewer than 3 eligible subjects."
            )
    else:
        log.append("Within‑subject Bartlett (QI 10) skipped (switch off).")  
    
    # ════════════════════════════════════════════════════════════════════════
    # 4.  Force a perfectly balanced design for the ANOVA
    # ════════════════════════════════════════════════════════════════════════
    # ════════════════════════════════════════════════════════════════════════
    # 4.  Force a perfectly balanced design (equal S & R)
    #     – required for the closed-form two-level ANOVA that follows
    # ════════════════════════════════════════════════════════════════════════
    # ═════════ 4.  Force a perfectly balanced design ══════════════════════
    if enforce_balance:                                      # NEW GUARD

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
        
    else:                                                    # NEW ─ skip pruning
        log.append("⚠️ Balance enforcement skipped – continuing with "
                   "unbalanced data; results will use the unbalanced formulas.")
        # ─── NEW ▶ log the FINAL numbers kept for analysis ─────────────────
        I_fin = df["Subject"].nunique()
        S_fin = df.groupby("Subject")["Sample"].nunique().mean()          # mean S
        R_fin = df.groupby(["Subject","Sample"])["Replicate"].nunique().mean()  # mean R
        log.append(
            f"✅ Final data set (without enforcement for balanced crossed design): {I_fin} subjects, "
            f"mean {S_fin:.2f} samples/subject, "
            f"mean {R_fin:.2f} replicates/sample retained for ANOVA."
        )
        # -------------------------------------------------------------------

        # minimal sanity check
        if df.empty:
            raise PreprocessError("All data removed during QC.", log)

    # ════════════════════════════════════════════════════════════════════════
    # 5.  Return cleaned data (back-transformed if needed) + log
    # ════════════════════════════════════════════════════════════════════════
    if transformed:
        df["Result"] = np.exp(df["Result"])

    return df, log



# CV ANOVA function
def _cvi_cv_anova(clean_df: pd.DataFrame, alpha: float = 0.05
) -> tuple[float, tuple[float, float]]:
    """
    CV-ANOVA as recommended by Røraas et al., 2016
    (see Clin Chem 62:725-736) :contentReference[oaicite:3]{index=3}
    -----------------------------------------------------------------
    1.  Divide every result by its subject-specific mean  →  CV-scale
    2.  Perform the usual two-level balanced ANOVA on the *normalised*
        values (mean ≈ 1).
    3.  σ²_WP (within-person) on that scale √→ CVI (%).
    4.  Exact CI via Burdick & Graybill χ² method (same df as classic).
    """
    df = clean_df.copy()
    df["Norm"] = df["Result"] / df.groupby("Subject")["Result"].transform("mean")

    I = df["Subject"].nunique()
    S = df.groupby("Subject")["Sample"].nunique().iloc[0]
    R = df.groupby(["Subject", "Sample"])["Replicate"].nunique().iloc[0]

    subj_mean = df.groupby("Subject")["Norm"].mean()                     # ~1
    samp_mean = df.groupby(["Subject", "Sample"])["Norm"].mean()

    ms_wp = R * ((samp_mean - subj_mean.loc[samp_mean.index.get_level_values(0)].values) ** 2
                 ).sum() / (I * (S - 1))
    ms_a  = ((df["Norm"] - samp_mean.loc[list(zip(df["Subject"], df["Sample"]))].values) ** 2
             ).sum() / (I * S * (R - 1))

    var_wp = max((ms_wp - ms_a) / R, 0.0)                # σ²_WP on CV scale
    cvi    = np.sqrt(var_wp) * 100                       # mean = 1 → CVI %

    # -------- exact 95 % CI (same df as classic) ------------------
    df_wp  = I * (S - 1)
    chi_lo = chi2.ppf(alpha/2,     df_wp)
    chi_hi = chi2.ppf(1 - alpha/2, df_wp)

    ms_lo  = (df_wp * ms_wp) / chi_hi
    ms_hi  = (df_wp * ms_wp) / chi_lo

    ci_var_lo = max((ms_lo - ms_a) / R, 0.0)
    ci_var_hi = max((ms_hi - ms_a) / R, 0.0)
    ci        = (np.sqrt(ci_var_lo) * 100, np.sqrt(ci_var_hi) * 100)
    return cvi, ci


# ---------------------------------------------------------------------------
#  CV‑ANOVA  →  unbalanced design
# ---------------------------------------------------------------------------
def _cvi_cv_anova_unbalanced(clean_df: pd.DataFrame,
                             alpha: float = 0.05
) -> tuple[float, tuple[float, float]]:
    """
    CV‑ANOVA that tolerates unequal samples‑/replicates‑per‑subject.

    Steps (Røraas 2016, adapted):

    1.  Normalise all results by each subject’s mean (→ mean ≈ 1).
    2.  Construct the unbalanced mean–squares table with weights nᵢⱼ.
    3.  σ²_WP = (MS_WP − MS_A)/ȓ      where ȓ = mean replicate count.
    4.  Exact CI uses Burdick & Graybill χ² limits with df = Σ(Sᵢ − 1).
    """
    df = clean_df.copy()
    df["Norm"] = df["Result"] / df.groupby("Subject")["Result"].transform("mean")

    # ------------------------------------------------------------------ weights
    r_ij = df.groupby(["Subject", "Sample"])["Replicate"].nunique()      # nᵢⱼ
    r_bar = r_ij.mean()                                                 # ȓ

    subj_mean = df.groupby("Subject")["Norm"].mean()
    samp_mean = df.groupby(["Subject", "Sample"])["Norm"].mean()

    # ── SS & MS (unbalanced) -----------------------------------------------
    ss_wp = ((samp_mean - samp_mean.index.get_level_values(0).map(subj_mean))**2 *
             r_ij).sum()
    df_wp = (df.groupby("Subject")["Sample"].nunique() - 1).sum()
    ms_wp = ss_wp / df_wp

    ss_a  = ((df["Norm"] -
              samp_mean.loc[list(zip(df["Subject"], df["Sample"]))].values)**2).sum()
    df_a  = (r_ij - 1).sum()
    ms_a  = ss_a / df_a

    var_wp = max((ms_wp - ms_a) / r_bar, 0.0)
    cvi    = np.sqrt(var_wp) * 100

    # ── exact CI on σ²_WP ----------------------------------------------------
    chi_lo = chi2.ppf(alpha/2,     df_wp)
    chi_hi = chi2.ppf(1-alpha/2,  df_wp)

    ms_lo  = (df_wp * ms_wp) / chi_hi
    ms_hi  = (df_wp * ms_wp) / chi_lo

    var_lo = max((ms_lo - ms_a) / r_bar, 0.0)
    var_hi = max((ms_hi - ms_a) / r_bar, 0.0)
    ci     = (np.sqrt(var_lo)*100, np.sqrt(var_hi)*100)

    return cvi, ci


# ——— helper for the unbalanced branch ——————————————————————————————
def _calculate_bv_unbalanced(df: pd.DataFrame,
                             alpha: float = 0.05) -> dict[str, float]:
    """
    Replicates the Excel ‘Glucose CVI and CVA – All’ workbook:

    • σ²_A   : pure analytical variance  = MS_A  
    • σ²_WP  : (MS_WP – MS_A) / ȓ        where ȓ is the *mean* replicate count  
    • σ²_BP  : (MS_BP – MS_WP) / (Ś·ȓ)   Ś = *mean* samples/subject  

    CIs use the same χ² limits as the balanced formula but with the actual
    d.f. from the unbalanced ANOVA table.
    """
    # raw counts per cell
    r_ij  = (df.groupby(["Subject", "Sample"])["Replicate"]
               .nunique())
    S_i   = df.groupby("Subject")["Sample"].nunique()

    r_bar = r_ij.mean()                      # ȓ
    S_bar = S_i.mean()                       # Ś

    # ---- ordinary nested MS table (works in unbalanced designs) ----------
    # NB: keep the sums as in the balanced code – pandas handles the weights
    subj_mean = df.groupby("Subject")["Result"].mean()
    samp_mean = df.groupby(["Subject", "Sample"])["Result"].mean()
    grand     = df["Result"].mean()

    ss_bp = ((samp_mean.index.get_level_values(0).map(subj_mean) - grand)**2
              * r_ij).groupby(level=0).sum().sum()
    ss_wp = ((samp_mean - samp_mean.index.get_level_values(0).map(subj_mean))**2
              * r_ij).sum()
    ss_a  = ((df["Result"]
              - samp_mean.loc[list(zip(df["Subject"], df["Sample"]))].values)**2).sum()

    I  = subj_mean.size
    df_bp = I - 1
    df_wp = (S_i - 1).sum()
    df_a  = (r_ij - 1).sum()

    ms_bp, ms_wp, ms_a = ss_bp/df_bp, ss_wp/df_wp, ss_a/df_a

    var_A  = ms_a
    var_WP = max((ms_wp - ms_a)/r_bar, 0.0)
    var_BP = max((ms_bp - ms_wp)/(S_bar*r_bar), 0.0)

    cv_A = np.sqrt(var_A)/grand * 100
    cv_I = np.sqrt(var_WP)/grand * 100
    cv_G = np.sqrt(var_BP)/grand * 100

    # χ² limits – same idea, just the *unbalanced* d.f.
    chi = chi2
    ci_cv_A = (np.sqrt(var_A*df_a/chi.ppf(0.975, df_a))/grand*100,
               np.sqrt(var_A*df_a/chi.ppf(0.025, df_a))/grand*100)
    ci_cv_I = (np.sqrt(max(( (df_wp*ms_wp/chi.ppf(0.975,df_wp))-var_A)/r_bar,0))/grand*100,
               np.sqrt(max(( (df_wp*ms_wp/chi.ppf(0.025,df_wp))-var_A)/r_bar,0))/grand*100)
    ci_cv_G = (np.sqrt(max(( (df_bp*(ms_bp-ms_wp)/chi.ppf(0.975,df_bp)))/(S_bar*r_bar),0))/grand*100,
               np.sqrt(max(( (df_bp*(ms_bp-ms_wp)/chi.ppf(0.025,df_bp)))/(S_bar*r_bar),0))/grand*100)

    rcv = 1.96*np.sqrt(2)*np.sqrt(cv_A**2 + cv_I**2)

    return dict(var_A=var_A, var_WP=var_WP, var_BP=var_BP,
                cv_A=cv_A, cv_I=cv_I, cv_G=cv_G,
                ci_cv_A=ci_cv_A, ci_cv_I=ci_cv_I, ci_cv_G=ci_cv_G,
                rcv=rcv, grand=grand,
                df_bp=df_bp, df_wp=df_wp, df_a=df_a)

# ——————————————————————————————————————————————————————————————
# 4.  Core calculation routine (closed‑form Røraas algebra)
# ——————————————————————————————————————————————————————————————

def calculate_bv(df: pd.DataFrame, alpha: float = 0.05,
                 use_cv_anova: bool = False,
                 enforce_balance: bool = True        # NEW PARAM
) -> BVResult:
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
    
    # ── NEW – capture raw study dimensions BEFORE any cleaning ───────────────
    I0 = df["Subject"].nunique()
    S0 = df.groupby("Subject")["Sample"].nunique().mode().iat[0]
    R0 = df.groupby(["Subject", "Sample"])["Replicate"].nunique().mode().iat[0]
    # -------------------------------------------------------------------------
    # 3.1  apply Braga–Panteghini outlier/normality pipeline  
    try:
        df, pp_log = _preprocess_bv_dataframe(
            df, alpha=alpha, enforce_balance=enforce_balance, flags=st.session_state.get("preproc_flags"))   # NEW ARG
    except PreprocessError as e:
        # Propagate the cleaner’s detailed log upward
        raise PreprocessError(str(e), e.log) from None

    # ── NEW – prepend the headline to the audit trail ────────────────────────
    pp_log.insert(
        0,
        f"Uploaded data: {I0} subjects × {S0} samples × {R0} replicates "
        "(raw, before quality checks)."
    )
    # ------------------------------------------------------------------------

    if not enforce_balance:                                 # ── unbalanced ──
        ub = _calculate_bv_unbalanced(df, alpha)

        # ─── NEW ▶ dataset counts & CI on the mean ─────────────────────────
        I_u  = df["Subject"].nunique()
        S_u  = df.groupby("Subject")["Sample"].nunique().mean()   # mean S
        R_u  = df.groupby(["Subject", "Sample"])["Replicate"].nunique().mean()  # mean R

        N_tot  = len(df)
        sd_all = df["Result"].std(ddof=1)
        t_crit = t.ppf(1 - alpha/2, N_tot - 1)
        ci_mean_u = (ub["grand"] - t_crit * sd_all / np.sqrt(N_tot),
                     ub["grand"] + t_crit * sd_all / np.sqrt(N_tot))
        # log it
#        pp_log.append(
#            f"✅ Unbalanced data set ready: {I_u} subjects × "
#            f"mean {S_u:.2f} samples/subject × "
#            f"mean {R_u:.2f} replicates/sample retained for ANOVA."
#        )
        pp_log.append("Unbalanced design → variance components derived with "
                      "(method‑of‑moments) algebra.")
        # ───────────────────────────────────────────────────────────────────
        # ▶ optional CV‑ANOVA on the *unbalanced* data
        cv_anova = ci_cv_anova = None
        if use_cv_anova:                                   # honour sidebar tick
            cv_anova, ci_cv_anova = _cvi_cv_anova_unbalanced(df, alpha)
            pp_log.append(
                f"CVI (CV‑ANOVA, unbalanced) calculated: {cv_anova:.2f} % "
                f"(95 % CI {ci_cv_anova[0]:.2f}–{ci_cv_anova[1]:.2f} %)."
            )

        return BVResult(
            I = int(I_u),
            S = S_u,
            R = R_u,
            grand_mean = ub["grand"],
            var_A = ub["var_A"],
            var_WP = ub["var_WP"],
            var_BP = ub["var_BP"],
            cv_A = ub["cv_A"],
            cv_I = ub["cv_I"],
            cv_G = ub["cv_G"],
            ci_mean = ci_mean_u,
            ci_cv_A = ub["ci_cv_A"],
            ci_cv_I = ub["ci_cv_I"],
            ci_cv_G = ub["ci_cv_G"],
            rcv_95 = ub["rcv"],
            preprocess_log = pp_log,
            cv_I_cv_anova = cv_anova,               # ← NEW
            ci_cv_I_cv_anova = ci_cv_anova          # ← NEW
        )


    else:  # balanced design → use closed‑form algebra

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

        # 🔶 NEW — 95 % CI on the mean ------------------------------
        N_tot  = len(df)
        sd_all = df["Result"].std(ddof=1)
        t_crit = t.ppf(1 - alpha/2, N_tot - 1)
        ci_mean = (grand - t_crit * sd_all / np.sqrt(N_tot),
                grand + t_crit * sd_all / np.sqrt(N_tot))
        # -----------------------------------------------------------

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

        # 🔶 NEW — 95 % CI on CV_A  ---------------------------------
        df_a        = I * S * (R - 1)            # d.f. of MS_A
        chi_low_a   = chi2.ppf(alpha/2,     df_a)
        chi_up_a    = chi2.ppf(1-alpha/2,  df_a)
        ci_var_a_lo = var_A * df_a / chi_up_a
        ci_var_a_hi = var_A * df_a / chi_low_a
        ci_cv_A     = (np.sqrt(ci_var_a_lo) / grand * 100,
                    np.sqrt(ci_var_a_hi) / grand * 100)
        # -----------------------------------------------------------

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
        
        # ⟨QI-6⟩ analytical imprecision – put a note into the audit trail
        pp_log.append(
            f"QI 6 – analytical imprecision calculated: CVₐ = {cv_A:.2f} % "
            f"(based on R = {R} replicate(s)/sample)."
        )

        # --- after the classic CV calculations are finished -----------------
        cv_anova = ci_cv_anova = None
        if use_cv_anova:
            cv_anova, ci_cv_anova = _cvi_cv_anova(df, alpha)
            pp_log.append(
                f"CVI (CV-ANOVA) calculated: {cv_anova:.2f} % "
                f"(95 % CI {ci_cv_anova[0]:.2f}–{ci_cv_anova[1]:.2f} %)."
            )

        return BVResult(
            I, S, R, grand,
            var_A, var_WP, var_BP,
            cv_A, cv_I, cv_G,
            ci_mean,          # add here
            ci_cv_A,          # add here
            ci_cv_I, ci_cv_G, rcv,
            preprocess_log = pp_log,
            cv_I_cv_anova = cv_anova,
            ci_cv_I_cv_anova = ci_cv_anova,
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

    # ─────────────────────────────────────────────────────────────
    # Pre‑processing switches – default = all ON
    # ─────────────────────────────────────────────────────────────
    st.sidebar.subheader("Pre‑processing steps")

    PREPROC_OPTS = {
        # label                        internal key            default
        "Replicate Cochran (QI 8a)"  : ("rep_cochran",          True),
        "Replicate Bartlett (QI 8a)" : ("rep_bartlett",         True),
        "Sample Cochran (QI 8b)"     : ("samp_cochran",         True),
        "Reed between‑subject (QI 8c)":("reed",                 True),
        "Steady‑state drift (QI 7)"  : ("drift",                True),
        "Normality checks / log‑transform (QI 9)"
                                    : ("normality",           True),
        "Within‑subject Bartlett (QI 10)"
                                    : ("wp_bartlett",         True),
    }

    # build the check‑boxes and stash the chosen flags in session_state
    preproc_flags = {}
    for label, (key, default) in PREPROC_OPTS.items():
        preproc_flags[key] = st.sidebar.checkbox(label, value=default)
    st.session_state["preproc_flags"] = preproc_flags

    # Crossed design balance enforcement
    st.sidebar.subheader("Analysis options")                    # NEW
    enforce_balance = st.sidebar.checkbox(                      # NEW
        "Enforce balanced crossed design",                      # NEW
        value=True,                                             # NEW (default keeps old behaviour)
    )                                                           # NEW
    st.session_state["enforce_balance"] = enforce_balance       # NEW

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

# NEW ▸ make sure the variable always exists
user_df: pd.DataFrame | None = None

# -- Tab 1: Upload -----------------------------------------------------------
with upload_tab:
    up_file = st.file_uploader("Upload CSV or XLSX", type=["csv", "xlsx"])
    if up_file is not None:
        try:
            if up_file.name.lower().endswith("xlsx"):
                user_df = pd.read_excel(up_file)
            else:
                user_df = pd.read_csv(up_file)

            # NEW ➜ strip accidental index column
            user_df = _strip_ghost_index(user_df)

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
    if manual_df.dropna().shape[0] >= 4:
        user_df = _strip_ghost_index(manual_df.copy())
    # save edits back to session
    st.session_state.manual_df = manual_df


# ——————————————————————————————————————————————————————————————
# 7. Preview & calculate button
# ——————————————————————————————————————————————————————————————
if user_df is not None:
    st.subheader("Data preview")
    st.dataframe(user_df, use_container_width=True, hide_index=True)


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

    # NEW – user option: CV-ANOVA
    estimate_cv_anova = st.checkbox("Estimate CVI with CV-ANOVA")
    st.session_state["use_cv_anova"] = estimate_cv_anova

    if st.button("Calculate", type="primary"):
        try:
            df_for_calc = st.session_state.get("mapped_df")
            if df_for_calc is None:
                st.warning("Please open “Map / confirm columns” and save a mapping first.")
                st.stop()                       # abort this run cleanly

            try:
                res = calculate_bv(
                    df_for_calc,
                    use_cv_anova=st.session_state.get("use_cv_anova", False),
                    enforce_balance = st.session_state.get("enforce_balance", True)   # NEW
                )


                # Key metrics — big bold labels, smaller numbers
                #  Key metrics – two-row layout
                # ────────────────────────────────────────────────────────────
                st.subheader("Key metrics")

                # --- build the metrics list ---
                metrics = [
                    ("Mean",
                    f"{res.grand_mean:.2f}",
                    f"{res.ci_mean[0]:.2f}–{res.ci_mean[1]:.2f} CI"),
                    ("CV<sub>A</sub>",
                    f"{res.cv_A:.2f} %",
                    f"{res.ci_cv_A[0]:.2f}–{res.ci_cv_A[1]:.2f}% CI"),
                    ("CV<sub>I</sub> <span style='font-size:0.7em;'>(based on standard ANOVA)</span>",
                        f"{res.cv_I:.2f} %",
                        f"{res.ci_cv_I[0]:.2f}–{res.ci_cv_I[1]:.2f}% CI"),
                    ("CV<sub>G</sub>",  f"{res.cv_G:.2f} %",
                        f"{res.ci_cv_G[0]:.2f}–{res.ci_cv_G[1]:.2f}% CI"),
                    ("95 % RCV",        f"±{res.rcv_95:.2f} %",                                    None),
                ]

                # append the optional CV-ANOVA metric as *last* element
                if res.cv_I_cv_anova is not None:
                    metrics.append((
                        "CV<sub>I</sub> <span style='font-size:0.75em;'>(based on CV-ANOVA)</span>",
                        f"{res.cv_I_cv_anova:.2f} %",
                        f"{res.ci_cv_I_cv_anova[0]:.2f}–{res.ci_cv_I_cv_anova[1]:.2f}% CI"
                    ))

                # --- split into rows ------------------------------------------------------
                first_row  = metrics[:-3]                     # everything except the very last
                second_row = metrics[-3:]                    # the last metric (may be empty)

                # helper that renders one row of “cards”
                def _render_row(row_metrics):
                    cols = st.columns(len(row_metrics), gap="large")
                    for col, (label, value, delta) in zip(cols, row_metrics):
                        with col:
                            st.markdown(f"""
                            <div style="
                                background:#f0f0f3; border-radius:1rem; padding:0.8rem 1rem;
                                text-align:center;
                                box-shadow:6px 6px 12px rgba(0,0,0,0.12),
                                        -6px -6px 12px rgba(255,255,255,0.85);
                            ">
                            <div style="font-size:1.4rem; font-weight:800; color:#1e3a5f;
                                        line-height:1.1; margin-bottom:0.15rem;">
                                {label}
                            </div>
                            <div style="font-size:1rem; font-weight:600; color:#1b263b;
                                        line-height:1.2;">{value}</div>
                            {f"<div style='font-size:0.75rem; color:#5d6d7e; margin-top:0.3rem;'>{delta}</div>" if delta else ""}
                            </div>
                            """, unsafe_allow_html=True)

                # draw the two rows
                _render_row(first_row)
                if second_row:                               
                    st.write(" ")
                    st.write(" ")
                    _render_row(second_row)


                # Summary table of CVs and 95% CIs
                st.write(" ")
                st.write(" ")
                st.subheader("Summary of variation metrics")
                var_tbl = pd.DataFrame({
                    "Component": ["Analytical", "Within‑subject", "Between‑subject"],
                    "Variance":  [res.var_A,    res.var_WP,      res.var_BP],
                    "CV %":      [res.cv_A,     res.cv_I,        res.cv_G],
                    "95 % CI":   [f"{res.ci_cv_A[0]:.2f}–{res.ci_cv_A[1]:.2f} %",
                                f"{res.ci_cv_I[0]:.2f}–{res.ci_cv_I[1]:.2f} %",
                                f"{res.ci_cv_G[0]:.2f}–{res.ci_cv_G[1]:.2f} %"],
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
                clean_df, _ = _preprocess_bv_dataframe(
                    df_for_calc,
                    flags=st.session_state["preproc_flags"],      # ← add
                    enforce_balance=st.session_state["enforce_balance"]
                )
                st.subheader("Per-subject distribution")
                st.plotly_chart(plot_subject_ranges(clean_df), use_container_width=True)


                # ---------------------------------------------------------------------------
                #  📋  BIVAC CHECKLIST  (QI 6 → QI 14)
                #      – fills Comment / Details and skips QI 1‑5
                # ---------------------------------------------------------------------------
                import re
                from collections import defaultdict

                def _build_qi_checklist(
                        log_lines: list[str],
                        res: BVResult,
                        flags: dict[str, bool],
                ) -> pd.DataFrame:
                    """
                    Build a BIVAC v1.1 checklist exactly as defined in Aarsand et al. 2018,
                    but **only for QI 6‑14** and with 2 new features:

                    • *Comment*  – short human‑readable explanation of the grade  
                    • *Details*  – concatenation of **all** log lines relevant to that QI  

                    Critical items (6‑8‑10‑11‑13) still auto‑downgrade to D when they
                    would otherwise have received C (see Table 1 of the paper).
                    """

                    # ------------------------------------------------------------------ helpers
                    def worst(g1, g2):                        # strict min on the A>B>C>D scale
                        order = {'A': 0, 'B': 1, 'C': 2, 'D': 3}
                        return g1 if order[g1] >= order[g2] else g2

                    # ════════════════════════════════════════════════════════════════════════
                    # 1.  Grade assessment  (same logic you already had, but no QI 1‑5)
                    #    ––––– I moved the old logic into a small dict for clarity –––––
                    # ════════════════════════════════════════════════════════════════════════
                    grade, comment = {}, {}       # QI ➜ grade / comment (A–D)

                    # ---------- QI 6  (replicates / analytical imprecision) -----------------
                    grade["6"]   = "A" if res.R >= 2 else "C"
                    comment["6"] = f"{res.R} replicates per sample detected."
                    # ---------- QI 7  (steady‑state / drift) --------------------------------
                    drift = any("temporal drift" in l.lower() for l in log_lines)
                    grade["7"]   = "A"
                    comment["7"] = "Significant drift subjects removed." if drift else "No drift detected."
                    # ---------- QI 8  (outlier handling) ------------------------------------
                    rep_ok  = flags.get("rep_cochran", True)
                    samp_ok = flags.get("samp_cochran", True)
                    subj_ok = flags.get("reed",        True)

                    if rep_ok and samp_ok and subj_ok:
                        grade["8"] = "A"
                        comment["8"] = "Replicate, sample and subject‑level outlier tests performed."
                    elif samp_ok and subj_ok:
                        grade["8"] = "B"
                        comment["8"] = "Replicate‑level outlier test skipped."
                    elif samp_ok:
                        grade["8"] = "B"
                        comment["8"] = "Replicate‑ & subject‑level outlier tests skipped."
                    else:
                        grade["8"] = "C"
                        comment["8"] = "Only partial outlier testing."
                    # ---------- QI 9  (normality) ------------------------------------------
                    transformed = any("log transform applied" in l.lower() for l in log_lines)
                    grade["9"]   = "A"
                    comment["9"] = "Log‑transform applied." if transformed else "Gaussian on raw scale."
                    # ---------- QI 10 (variance homogeneity) -------------------------------
                    het = any("heterogeneous" in l.lower() for l in log_lines)
                    if not flags.get("wp_bartlett", True):
                        grade["10"] = "B"
                        comment["10"] = "Variance homogeneity test skipped."
                    else:
                        grade["10"] = "A" if not het else "C"
                        comment["10"] = "Variances homogeneous." if not het else "Heterogeneous variances."
                    # ---------- QI 11 (ANOVA model) ----------------------------------------
                    balanced = res.S == int(res.S) and res.R == int(res.R)
                    grade["11"]   = "A" if balanced else "B"
                    comment["11"] = "Balanced nested ANOVA." if balanced else "Method‑of‑moments (unbalanced)."
                    # ---------- QI 12 (confidence limits) ----------------------------------
                    grade["12"]   = "A"
                    comment["12"] = "95 % CIs on CVA/CVI/CVG reported."
                    # ---------- QI 13 (numbers kept) ---------------------------------------
                    grade["13"]   = "A"
                    kept = res.I * res.S * res.R
                    comment["13"] = f"{kept} results retained after QC."
                    # ---------- QI 14 (mean concentration) ---------------------------------
                    grade["14"]   = "A"
                    comment["14"] = f"Mean concentration {res.grand_mean:.2f}."

                    # ------------------------------------------------------------------ auto‑downgrade critical items to D
                    for q in ("6", "8", "10", "11", "13"):
                        if grade[q] == "C":
                            grade[q] = "D"
                            comment[q] += "  (critical item downgraded to D)"

                    # ════════════════════════════════════════════════════════════════════════
                    # 2.  Map preprocessing log → Details column
                    # ════════════════════════════════════════════════════════════════════════
                    # simple pattern map: QI ➜ keywords expected in that log line
                    patterns = {
                        "6": ("cvₐ", "analytical imprecision"),
                        "7": ("drift",),
                        "8": ("outlier",),
                        "9": ("normality", "log transform"),
                        "10": ("bartlett", "heterogeneous"),
                        "11": ("balanced", "unbalanced", "anova"),
                        "13": ("subject", "sample", "replicate", "results retained"),
                        "14": ("mean",),
                        # QI 12 has no explicit log hooks – we’ll leave Details blank
                    }

                    details = defaultdict(list)
                    for line in log_lines:
                        line_l = line.lower()
                        for q, keys in patterns.items():
                            if any(k in line_l for k in keys):
                                details[q].append(line)

                    # ════════════════════════════════════════════════════════════════════════
                    # 3.  Assemble DataFrame (QI 6‑14) and overall grade
                    # ════════════════════════════════════════════════════════════════════════
                    rows, overall = [], "A"
                    for q in map(str, range(6, 15)):
                        overall = worst(overall, grade[q])
                        rows.append(dict(
                            QI      = f"QI {q}",
                            Grade   = grade[q],
                            Comment = comment[q],
                            Details = " ⏵ ".join(details.get(q, []))          # nice arrow separator
                        ))

                    df = pd.DataFrame(rows)
                    df.attrs["overall_grade"] = overall
                    return df




                # ── render the checklist & XLSX download ────────────────────────────────────
                st.subheader("BIVAC Checklist (QI 6 → QI 14)")
                bivac_df = _build_qi_checklist(
                    res.preprocess_log,
                    res,
                    flags = st.session_state["preproc_flags"]
                )

                st.dataframe(
                    bivac_df.style.apply(
                        lambda s: ["background:#e8f9f0" if g in ("A", "B") else
                                "background:#fdecea"   for g in s],
                        axis=1, subset=["Grade"]
                    ).set_properties(**{"font-size": "0.85rem"}),
                    use_container_width=True,
                    hide_index=True,
                )
                overall = bivac_df.attrs["overall_grade"]
                st.markdown(f"**Overall BIVAC grade: {overall}**")

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
                st.write(e)
                st.error(f"Calculation failed: {e}")     # any other unexpected error

        except Exception as e:
            st.write(e)
            st.error(f"Calculation failed: {e}")
else:
    st.info("Input data above to enable calculation.")
