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
from scipy.stats import shapiro, f              # normality + Cochran
from scipy.stats import bartlett, linregress  
import plotly.graph_objects as go  
import re, html


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
      /* ===== Key metrics (table) – fancy HTML table ===== */
      .km-wrap { margin-top: .5rem; }
      table.km-table {
        width: 100%;
        border-collapse: separate;
        border-spacing: 0 .45rem;   /* “card rows” */
        font-family: 'Georgia', serif;
      }
      table.km-table thead th{
        text-align: left;
        font-size: .85rem;
        letter-spacing: .03em;
        color: #1e3a5f;
        background: #d0e7f5;
        padding: .65rem .8rem;
        border: 1px solid #b0c4de;
      }
      table.km-table thead th:first-child{ border-top-left-radius: .8rem; border-bottom-left-radius: .8rem; }
      table.km-table thead th:last-child { border-top-right-radius:.8rem; border-bottom-right-radius:.8rem; }

      table.km-table tbody tr{
        background: linear-gradient(145deg, #ffffff, #f6fbff);
        box-shadow: 0 2px 8px rgba(0,0,0,0.05);
      }
      table.km-table tbody td{
        padding: .75rem .8rem;
        border-top: 1px solid #e3eef8;
        border-bottom: 1px solid #e3eef8;
      }
      table.km-table tbody td:first-child{
        border-left: 1px solid #e3eef8;
        border-top-left-radius: .85rem;
        border-bottom-left-radius: .85rem;
        font-weight: 700;
        color: #1e3a5f;
      }
      table.km-table tbody td:last-child{
        border-right: 1px solid #e3eef8;
        border-top-right-radius: .85rem;
        border-bottom-right-radius: .85rem;
        text-align: right;
        white-space: nowrap;
      }
      table.km-table tbody tr:hover{
        transform: translateY(-1px);
        box-shadow: 0 6px 14px rgba(0,0,0,0.08);
      }

      .km-pill{
        display:inline-block;
        padding: .12rem .5rem;
        border-radius: 999px;
        background: #e8f4fd;
        border: 1px solid #a3cce3;
        color: #28527a;
        font-size: .72rem;
        margin-right: .45rem;
        vertical-align: middle;
      }
      .km-ci{
        color:#5d6d7e;
        font-size:.78rem;
        margin-left:.4rem;
        white-space: nowrap;
      }
      .km-value{
        font-weight: 800;
        color:#1b263b;
      }
      /* ==== Preprocessing log highlighting ==== */
    .log-list{ list-style: none; padding-left: 0; margin: 0; }
    .log-item{ margin: .25rem 0; }
    .log-outlier{ color:#b00020; font-weight:700; }   /* red = outliers/exclusions */
    .log-sig{     color:#a15c00; font-weight:700; }   /* amber = significant stats */
    .log-action{  color:#004c97; font-weight:700; }   /* blue = major manipulation */
    </style>
    """,
    unsafe_allow_html=True
)
# ---- Always define these so later guards don't crash ----
user_df: pd.DataFrame | None = None
up_file = None

# --- BIVAC QI 1–5 definitions (Aarsand 2018 Table 1) ---
QI15_CHOICES = {
    "QI 1 – Scale": {
        "help": " Is the analyte expressed on a ratio scale?",
        "options": {
            "A": "Yes (ratio scale)",
            "B": "No (non-ratio scale)",
        },
    },
    "QI 2 – Subjects": {
        "help": "Subjects/population documented: (a) number; (b) sex; (c) age/age group; (d) health status.",
        "options": {
            "A": "All (a,b,c,d) documented",
            "B": "All documented but only 'healthy volunteers' given",
            "C": "(a) and (d) documented; missing (b)/(c) not important",
            "D": "(a) and (d) documented; missing (b)/(c) important OR (a)/(d) missing",
        },
    },
    "QI 3 – Samples": {
        "help": "Are these aspects recorded? (a) Count of collected samples; (b) Specimen type used; (c) Collection timing; (d) Study duration",
        "options": {
            "A": "Yes (a,b,c,d) documented",
            "B": "(a,c,d) documented; (b) insufficient but not important",
            "C": "(a,c,d) documented; (b) insufficient and important",
            "D": "(a), (c) and/or (d) not presented/deducible",
        },
    },
    "QI 4 – Measurand/Method": {
        "help": "Method documented sufficiently / fit for BV estimation.",
        "options": {
            "A": "Detailed method or adequate reference/identifiable method",
            "B": "Insufficient detail (not important for measurand)",
            "C": "Insufficient detail or outdated (important for measurand)",
            "D": "Obsolete / not fit for BV estimation",
        },
    },
    "QI 5 – Preanalytical": {
        "help": "Preanalytical procedures standardized to minimize variation.",
        "options": {
            "A": "Yes (adequate)",
            "B": "Insufficient detail (unlikely important)",
            "C": "Insufficient detail (may be important) or No details provided",
        },
    },
    "QI 6 – Analytical Variation": {
        "help": "Estimates of analytical variation: same-run vs. different-run replicates.",
        "options": {
            "A": "All replicates for the same subject analyzed in the same run",
            "B": "Replicate analyses of samples performed in different runs",
        },
    },
}



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
    rcv_95_down: float   # % decrease (reported as a positive number)
    rcv_95_up: float     # % increase

    # NEW – verbose audit trail of what happened during cleaning
    preprocess_log: list[str] = field(default_factory=list)
    
    # NEW – optional CV-ANOVA estimates
    cv_I_cv_anova: float | None = None
    ci_cv_I_cv_anova: tuple[float, float] | None = None

    # NEW: store the exact cleaned dataframe used for the final ANOVA
    clean_df: pd.DataFrame | None = None
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

# ─────────────────────────────────────────────────────────────────────────────
# Mapping validation
# ─────────────────────────────────────────────────────────────────────────────
def _validate_and_build_mapping(df: pd.DataFrame,
                                subj_sel: str, samp_sel: str, repl_sel: str, res_sel: str
) -> tuple[pd.DataFrame | None, list[str], list[str], list]:
    """Return (mapped_df, errors, warnings, dropped_indices). Ensures 4 distinct columns and numeric 'Result'."""
    errors, warnings = [], []
    dropped_indices = []  # Track indices of dropped rows for gender-specific reporting

    # Distinct columns?
    chosen = [subj_sel, samp_sel, repl_sel, res_sel]
    if len(set(chosen)) != 4:
        errors.append("Please select **four different columns** (no duplicates).")

    # Build mapped view early so we can check types - preserve original index
    mapped = df.rename(columns={
        subj_sel: "Subject",
        samp_sel: "Sample",
        repl_sel: "Replicate",
        res_sel:  "Result",
    })[["Subject", "Sample", "Replicate", "Result"]].copy()

    # Check for missing/empty Result values BEFORE numeric conversion
    # Count rows where Result is NaN, None, or empty string
    original_result = mapped["Result"]
    missing_mask = original_result.isna() | (original_result.astype(str).str.strip() == "")
    n_missing = missing_mask.sum()

    # Coerce numeric columns
    mapped["Result"] = pd.to_numeric(mapped["Result"], errors="coerce")
    n_na_after = mapped["Result"].isna().sum()

    # Non-numeric = values that became NaN after conversion but weren't originally missing
    n_non_numeric = n_na_after - n_missing

    if n_non_numeric > 0:
        errors.append(f"'Result / value' must be numeric → {n_non_numeric} row(s) contain non-numeric text.")

    # Handle missing values: drop them with a warning instead of error
    if n_missing > 0:
        warnings.append(f"Dropped {n_missing} row(s) with missing/empty Result values.")
        # Store indices of dropped rows before dropping
        dropped_indices = mapped[mapped["Result"].isna()].index.tolist()
        mapped = mapped.dropna(subset=["Result"])

    # Sample & Replicate should be numeric (visit/timepoint & within-sample repeat index)
    mapped["Sample"] = pd.to_numeric(mapped["Sample"], errors="coerce")
    mapped["Replicate"] = pd.to_numeric(mapped["Replicate"], errors="coerce")
    if mapped["Sample"].isna().any():
        errors.append("'Sample / time-point' must be numeric (e.g., 1, 2, 3).")
    if mapped["Replicate"].isna().any():
        errors.append("'Replicate' must be numeric (e.g., 1, 2).")

    # Light sanity notes
    if mapped["Subject"].isna().any():
        warnings.append("Some Subject IDs are blank.")
    if mapped.dropna().empty:
        errors.append("All mapped values are empty after type checks.")

    return (None, errors, warnings, dropped_indices) if errors else (mapped, errors, warnings, dropped_indices)


def _render_log_html(lines: list[str], alpha: float = 0.05) -> str:
    """Return an HTML bullet list where important log messages are bold & colored."""
    def classify(line: str) -> str:
        L = line.lower()

        # --- Make "no outlier ..." lines plain (early exit) ---
        if ("no outlier detected" in L
            or "no between-subject outliers" in L
            or "no outliers detected" in L):
            return ""

        # --- Heterogeneous anywhere → emphasize (amber) ---
        if "heterogeneous" in L:
            return "log-sig"   # bold + amber via your CSS

        # --- Outliers / removals (red) ---
        if ("outlier removed" in L or "reed mean outlier" in L or " → outlier" in L
            or ("excluded" in L and "drift" in L)):
            return "log-outlier"

        # --- Significant tests (amber) ---
        m = re.search(r"\bp\s*=\s*([0-9]*\.?[0-9]+)", L)
        if m:
            try:
                if float(m.group(1)) < alpha:
                    return "log-sig"
            except ValueError:
                pass
        if ("gcrit" in L and "→ outlier" in L):
            return "log-sig"

        if "high within-subject variance" in L or "high analytical variance" in L:
            return "log-sig"

        # --- Major manipulations (blue) ---
        if any(k in L for k in [
            "natural-log transform applied", "log transform applied",
            "balance check", "subject **excluded**", "sample **removed**",
            "dropped", "discarded",
            "unbalanced design →", "balanced data set ready",
        ]):
            return "log-action"

        # otherwise: plain
        return ""

    items = []
    for ln in lines:
        cls = classify(ln)
        safe = html.escape(ln)
        if cls:
            items.append(f"<li class='log-item'><span class='{cls}'>{safe}</span></li>")
        else:
            items.append(f"<li class='log-item'>{safe}</li>")
    return "<ul class='log-list'>" + "\n".join(items) + "</ul>"

def _style_details_series(s: pd.Series) -> list[str]:
    """Make 'Details' text non-colored; keep bold for important rows."""
    triggers = [
        "outlier", "excluded", "dropped", "discarded",
        "log transform", "heterogeneous", "temporal drift",
        "high within-subject variance", "high analytical variance"
    ]
    styles = []
    for val in s.astype(str).str.lower():
        if any(k in val for k in triggers):
            styles.append("font-weight:700;")          # no color
        else:
            styles.append("")                           # default text color
    return styles

# population trend plot
def plot_population_trend(clean_df: pd.DataFrame) -> go.Figure:
    """
    Show cohort mean at each Sample (time) with an OLS fit line.
    X-axis ticks are forced to show *every* Sample value.
    """
    # collapse to Subject×Sample means, then mean across subjects per Sample
    ds = (clean_df.groupby(["Subject", "Sample"], as_index=False)["Result"].mean())
    ts = (ds.groupby("Sample")["Result"].agg(mean="mean", count="count", std="std").reset_index())

    # make sure sample axis is numeric and sorted
    samples = ts["Sample"].astype(int).sort_values().to_numpy()
    means   = ts.set_index("Sample").loc[samples, "mean"].to_numpy()

    # simple OLS fit on (Sample, mean)
    if len(samples) >= 2:
        b1, b0 = np.polyfit(samples.astype(float), means.astype(float), 1)  # slope, intercept
        yhat = b1 * samples + b0
    else:
        yhat = means  # not enough points to fit

    unit = st.session_state.get("result_unit", "").strip() if "result_unit" in st.session_state else ""
    ylab = "Result" + (f" ({unit})" if unit else "")

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=samples, y=means, mode="markers", name="Mean per time"))
    fig.add_trace(go.Scatter(x=samples, y=yhat,   mode="lines",   name="OLS fit"))

    # Force a tick at *every* sample value
    tickvals = samples.tolist()
    fig.update_layout(
        template="simple_white",
        height=280,
        margin=dict(l=40, r=10, t=30, b=40),
        xaxis=dict(
            title="Sample (time)",
            tickmode="array",
            tickvals=tickvals,
            ticktext=[str(v) for v in tickvals],
            range=[samples.min() - 0.5, samples.max() + 0.5],  # small padding left/right
            showgrid=False
        ),
        yaxis=dict(title=ylab, zeroline=False, showgrid=False),
        showlegend=False
    )
    return fig


def _build_qi_checklist(
        log_lines: list[str],
        res_: BVResult,
        flags: dict[str, bool],
        *,
        qi_manual: dict[str, str] | None = None,
        n_raw: int | None = None,
        n_kept: int | None = None,
        unit: str = "",
) -> pd.DataFrame:
    """BIVAC v1.1 checklist (QI 1–14) with Grade, Comment, and Details."""

    def worst(g1, g2):
        order = {'A': 0, 'B': 1, 'C': 2, 'D': 3}
        return g1 if order[g1] >= order[g2] else g2

    def pick(*keywords):
        hits = [l for l in log_lines if all(k.lower() in l.lower() for k in keywords)]
        return hits

    def any_kw(*keywords):
        return any(k.lower() in l.lower() for l in log_lines for k in keywords)

    details = {str(i): [] for i in range(1, 15)}
    grade, comment = {}, {}

    # ---------------- QI 1–5 (manual) ----------------
    if qi_manual:
        g1 = qi_manual.get("QI 1 – Scale", "A")
        g2 = qi_manual.get("QI 2 – Subjects", "A")
        g3 = qi_manual.get("QI 3 – Samples", "A")
        g4 = qi_manual.get("QI 4 – Measurand/Method", "A")
        g5 = qi_manual.get("QI 5 – Preanalytical", "A")
        grade.update({"1": g1, "2": g2, "3": g3, "4": g4, "5": g5})
        comment.update({
            "1": f"User-selected {g1}.",
            "2": f"User-selected {g2}.",
            "3": f"User-selected {g3}.",
            "4": f"User-selected {g4}.",
            "5": f"User-selected {g5}."
        })
        details["1"].append(f"Manual: Scale → {g1}.")
        details["2"].append(f"Manual: Subjects → {g2}.")
        details["3"].append(f"Manual: Samples → {g3}.")
        details["4"].append(f"Manual: Measurand/Method → {g4}.")
        details["5"].append(f"Manual: Preanalytical → {g5}.")

    # ---------------- QI 6 (Estimates of analytical variation) ----------------
    # BIVAC Table 1:
    # A: Estimates presented, all replicates for same subject analyzed in same run
    # B: Estimates presented but obtained by other method OR replicates in different runs
    # C: No estimates presented
    qi6_grade = (qi_manual or {}).get("QI 6 – Analytical Variation", "A")
    if res_.R >= 2:
        grade["6"] = qi6_grade  # User selection (A=same run, B=different runs)
        if qi6_grade == "A":
            comment["6"] = f"Estimates presented ({res_.R} replicates/sample); same run analysis."
        else:
            comment["6"] = f"Estimates presented ({res_.R} replicates/sample); different runs."
    else:
        grade["6"] = "C"
        comment["6"] = f"No replicate analysis (only {res_.R} replicate/sample)."
    details["6"].append(f"User selected: Grade {qi6_grade}.")
    details["6"].append(f"CVₐ {res_.cv_A:.2f}% (95% CI {res_.ci_cv_A[0]:.2f}–{res_.ci_cv_A[1]:.2f}%).")

    # ---------------- QI 7 (Steady state) ----------------
    # BIVAC Table 1:
    # A: Yes - individual trend analysis performed or data adequately transformed
    # B: Not performed, but unlikely important for measurand
    # C: Not performed, and may be important (e.g., hormones) or clinical setting
    drift_checked = flags.get("drift", True)
    drift_removed = any_kw("temporal drift", "removed") or any_kw("drift outlier removed")
    drift_passed = any_kw("Steady-state drift test passed")
    
    if drift_checked:
        grade["7"] = "A"
        if drift_removed:
            comment["7"] = "Individual trend analysis performed; drifting subjects removed."
        elif drift_passed:
            comment["7"] = "Individual trend analysis performed; all subjects in steady state."
        else:
            comment["7"] = "Individual trend analysis performed."
    else:
        # Use user selection from sidebar (qi7_importance)
        qi7_grade = flags.get("qi7_importance", "C")  # Default to C if not set
        if qi7_grade not in ["A", "B", "C"]:
            qi7_grade = "C"
        grade["7"] = qi7_grade
        if qi7_grade == "B":
            comment["7"] = "Individual trend analysis not performed; unlikely important for this measurand."
        else:
            comment["7"] = "Individual trend analysis not performed; may be important for this measurand."
    details["7"].extend(pick("QI 7", "drift"))
    details["7"].extend(pick("temporal drift"))
    if not details["7"]:
        details["7"].append("No drift messages in log.")

    # ---------------- QI 8 (Outliers) ----------------
    # BIVAC Table 1:
    # A: Testing for outliers of (a) Replicates, (b) Samples per subject, (c) Subjects - ALL performed
    # B: (b) Fulfilled, but (a) replicate analysis not performed, and/or (c) subject outlier analysis not performed
    # C: (b) Sample-level outlier analysis not performed, or only performed on total data set
    rep_ok = flags.get("rep_cochran", True)
    samp_ok = flags.get("samp_cochran", True)
    subj_ok = flags.get("reed", True)
    
    if rep_ok and samp_ok and subj_ok:
        grade["8"] = "A"
        comment["8"] = "Outlier testing performed: (a) Replicates, (b) Samples, (c) Subjects."
    elif samp_ok:
        grade["8"] = "B"
        missing = []
        if not rep_ok:
            missing.append("replicate")
        if not subj_ok:
            missing.append("subject")
        comment["8"] = f"Sample-level outlier testing performed; {'/'.join(missing)} outlier testing not performed."
    else:
        grade["8"] = "C"
        comment["8"] = "Sample-level outlier analysis not performed."
    
    details["8"].extend(pick("Cochran", "replicate"))
    details["8"].extend(pick("Cochran", "sample"))
    details["8"].extend(pick("Reed"))
    if not details["8"]:
        details["8"].append("No outlier messages in log.")

    # ---------------- QI 9 (Normally distributed data) ----------------
    # BIVAC Table 1:
    # A: Yes - distribution assessed for each subject, transformed if not normal
    # B: No - distribution not assessed
    norm_checked = flags.get("normality", True)
    log_transformed = any_kw("log transform")
    
    if norm_checked:
        grade["9"] = "A"
        if log_transformed:
            comment["9"] = "Distribution assessed; log-transformation applied."
        else:
            comment["9"] = "Distribution assessed for normality."
    else:
        grade["9"] = "B"
        comment["9"] = "Distribution not assessed."
    
    details["9"].extend(pick("Normality"))
    details["9"].extend(pick("Shapiro-Wilk"))
    details["9"].extend(pick("log transform"))
    if not details["9"]:
        details["9"].append("No normality messages in log.")

    # ---------------- QI 10 (Variance homogeneity) ----------------
    # BIVAC Table 1:
    # A: Yes - variance homogeneity examined
    # C: No - variance homogeneity not examined
    wp_examined = flags.get("wp_bartlett", True)
    rep_bartlett_examined = flags.get("rep_bartlett", True)
    variance_homogeneity_examined = wp_examined or rep_bartlett_examined
    
    if variance_homogeneity_examined:
        grade["10"] = "A"
        het = any(("heterogeneous" in l.lower()) and ("bartlett" in l.lower()) for l in log_lines)
        hom = any(("homogeneous" in l.lower()) and ("heterogeneous" not in l.lower()) and ("bartlett" in l.lower()) for l in log_lines)
        if het:
            comment["10"] = "Variance homogeneity examined; heterogeneous subjects handled."
        elif hom:
            comment["10"] = "Variance homogeneity examined; variances are homogeneous."
        else:
            comment["10"] = "Variance homogeneity examined."
    else:
        grade["10"] = "C"
        comment["10"] = "Variance homogeneity not examined."
    
    details["10"].extend(pick("Bartlett"))
    details["10"].extend(pick("variance homogeneity"))
    if not details["10"]:
        details["10"].append("No variance-homogeneity messages in log.")

    # ---------------- QI 11 (Statistical method) ----------------
    # BIVAC Table 1:
    # A: Nested ANOVA or equivalent variance decomposition with estimation of analytical, within-, and between-subject variation
    # B: Simple subtraction of variances used
    # C: Other method or method not declared
    grade["11"] = "A"
    if (res_.S == int(res_.S)) and (res_.R == int(res_.R)):
        comment["11"] = "Nested ANOVA with variance decomposition (balanced design)."
        details["11"].append("Balanced crossed design; closed-form variance decomposition for CVₐ, CVᵢ, CVg.")
    else:
        comment["11"] = "Variance decomposition (unbalanced data; method-of-moments estimation)."
        details["11"].append("Unbalanced design; method-of-moments variance decomposition for CVₐ, CVᵢ, CVg.")

    # ---------------- QI 12 (Confidence limits) ----------------
    # BIVAC Table 1:
    # A: Yes - confidence limits around estimates of CVI presented
    # C: No - confidence limits not presented
    ci_available = (res_.ci_cv_I is not None and len(res_.ci_cv_I) == 2)
    if ci_available:
        grade["12"] = "A"
        comment["12"] = "95% confidence limits for CVI presented."
        details["12"].append(
            f"CIs: CVₐ {res_.ci_cv_A[0]:.2f}–{res_.ci_cv_A[1]:.2f}%, "
            f"CVᵢ {res_.ci_cv_I[0]:.2f}–{res_.ci_cv_I[1]:.2f}%, "
            f"CVg {res_.ci_cv_G[0]:.2f}–{res_.ci_cv_G[1]:.2f}%."
        )
    else:
        grade["12"] = "C"
        comment["12"] = "Confidence limits not available."
        details["12"].append("CIs could not be calculated.")

    # ---------------- QI 13 (Number of results) ----------------
    # BIVAC Table 1:
    # A: Number of results given and exclusion criteria defined with number of exclusions
    # B: Number of results and exclusion criteria given, but number excluded not stated
    # C: Number of results not given, or exclusion criteria not given
    if (n_raw is not None) and (n_kept is not None):
        n_excl = max(n_raw - n_kept, 0)
        grade["13"] = "A"
        comment["13"] = f"Results: {n_kept} used / {n_raw} raw ({n_excl} excluded)."
        details["13"].append(f"Kept {n_kept} of {n_raw} measurements; {n_excl} excluded by QC.")
    elif n_kept is not None:
        grade["13"] = "B"
        comment["13"] = f"{n_kept} results used; exclusion count not documented."
    else:
        grade["13"] = "C"
        comment["13"] = "Number of results not documented."

    # ---------------- QI 14 (Mean/median concentration) ----------------
    # BIVAC Table 1:
    # A: Yes - mean/median concentration reported with reference interval
    # B: Mean/median reported, reference interval not given
    # C: Mean/median not reported
    if res_.grand_mean is not None:
        grade["14"] = "A"  # Always A if we have the mean (reference interval from CIs)
        comment["14"] = (
            f"Mean concentration: {res_.grand_mean:.2f}"
            + (f" {unit}" if unit else "")
            + f" (95% CI {res_.ci_mean[0]:.2f}–{res_.ci_mean[1]:.2f}"
            + (f" {unit}" if unit else "")
            + ")."
        )
        details["14"].append(
            f"Mean = {res_.grand_mean:.2f}"
            + (f" {unit}" if unit else "")
            + f"; 95% CI = {res_.ci_mean[0]:.2f}–{res_.ci_mean[1]:.2f}"
            + (f" {unit}" if unit else "")
            + "."
        )
    else:
        grade["14"] = "C"
        comment["14"] = "Mean concentration not available."

    rows, overall = [], "A"
    for q in map(str, range(1, 15)):
        overall = worst(overall, grade[q])
        rows.append(dict(
            QI=f"QI {q}",
            Grade=grade[q],
            Comment=comment[q],
            Details=" ⏵ ".join(details[q]) if details[q] else ""
        ))
    out = pd.DataFrame(rows)
    out.attrs["overall_grade"] = overall
    return out





def build_bivac_df_for_group(label: str, r: BVResult, raw_df: pd.DataFrame) -> pd.DataFrame:
    final_df_x = r.clean_df
    if final_df_x is None or final_df_x.empty:
        raise ValueError(f"{label}: cleaned dataset is empty after QC.")

    unit = st.session_state.get("result_unit", "").strip()
    flags = st.session_state.get("preproc_flags", {})
    qi_manual = st.session_state.get("qi_manual", None)

    bivac_df = _build_qi_checklist(
        r.preprocess_log,
        r,
        flags=flags,
        qi_manual=qi_manual,
        n_raw=len(raw_df),
        n_kept=len(final_df_x),
        unit=unit,
    )
    return bivac_df


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


def _detect_replicate_outliers(df: pd.DataFrame, alpha: float = 0.05) -> list[dict]:
    """
    Detect replicate-level Cochran outliers WITHOUT removing them.
    Returns a list of dicts with outlier information for user selection.
    
    Each dict contains:
        - Subject: subject ID
        - Sample: sample ID  
        - replicates: list of {Replicate, Result} dicts
        - variance: the sample variance of the replicate set
        - G: Cochran G statistic
        - G_crit: Cochran critical value
    """
    outliers = []
    df_work = df.copy()
    
    while True:
        rep_var = (df_work.groupby(["Subject", "Sample"])["Result"]
                   .var(ddof=1).dropna())
        
        if rep_var.size < 3:
            break
        
        R = (df_work.groupby(["Subject", "Sample"])["Replicate"]
             .nunique().mode().iat[0])
        
        if R < 2:
            break
        
        c_res = _cochrans_test(rep_var.values, alpha, df=R - 1)
        
        if not c_res.flag:
            break
        
        subj, samp = rep_var.idxmax()
        rep_data = df_work.loc[(df_work["Subject"] == subj) & (df_work["Sample"] == samp),
                               ["Replicate", "Result"]].to_dict("records")
        
        outliers.append({
            "Subject": subj,
            "Sample": samp,
            "replicates": rep_data,
            "variance": float(rep_var.max()),
            "G": float(c_res.G),
            "G_crit": float(c_res.G_crit),
        })
        
        # Remove this pair from working copy to find next outlier
        df_work = df_work[~((df_work["Subject"] == subj) & (df_work["Sample"] == samp))]
    
    return outliers


def _detect_bartlett_outliers(df: pd.DataFrame, alpha: float = 0.05) -> list[dict]:
    """
    Detect Bartlett-based replicate outliers WITHOUT removing them.
    Returns a list of dicts with outlier information for user selection.
    
    Similar to Cochran detection, but uses the Bartlett test.
    Each dict contains:
        - Subject: subject ID
        - Sample: sample ID  
        - replicates: list of {Replicate, Result} dicts
        - variance: the sample variance of the replicate set
        - chi2: Bartlett chi-squared statistic
        - chi2_crit: critical value
        - p: p-value
    """
    outliers = []
    df_work = df.copy()
    
    while True:
        valid_groups = []
        group_keys = []
        for (subj, samp), g in df_work.groupby(["Subject", "Sample"]):
            if g["Result"].size < 2:
                continue
            if g["Result"].var(ddof=1) == 0:
                continue
            valid_groups.append(g["Result"].values)
            group_keys.append((subj, samp))
        
        if len(valid_groups) < 3:
            break
        
        with np.errstate(all="ignore"):
            bart_stat, bart_p = bartlett(*valid_groups)
        
        k = len(valid_groups)
        chi_crit = chi2.ppf(1 - alpha, k - 1)
        
        if bart_p >= alpha:
            # Homogeneous - done
            break
        
        # Heterogeneous - find worst offender
        rep_var = (df_work.groupby(["Subject", "Sample"])["Result"]
                   .var(ddof=1)
                   .loc[group_keys])
        
        (subj, samp) = rep_var.idxmax()
        max_var = rep_var.max()
        
        rep_data = df_work.loc[(df_work["Subject"] == subj) & (df_work["Sample"] == samp),
                               ["Replicate", "Result"]].to_dict("records")
        
        outliers.append({
            "Subject": subj,
            "Sample": samp,
            "replicates": rep_data,
            "variance": float(max_var),
            "chi2": float(bart_stat),
            "chi2_crit": float(chi_crit),
            "p": float(bart_p),
        })
        
        # Remove this pair from working copy to find next outlier
        df_work = df_work[~((df_work["Subject"] == subj) & (df_work["Sample"] == samp))]
    
    return outliers


def _detect_sample_cochran_outliers(df: pd.DataFrame, alpha: float = 0.05) -> list[dict]:
    """
    Detect Sample-level Cochran outliers (QI 8b) within each subject.
    Returns a list of dicts with outlier information for user selection.
    
    Iteratively detects outliers but does NOT permanent remove them from input df.
    """
    outliers = []
    
    # Process each subject separately
    for subj, g_subj in df.groupby("Subject"):
        g_work = g_subj.copy()
        
        while True:
            # Calculate variance for each sample within subject
            samp_var = g_work.groupby("Sample")["Result"].var(ddof=1).dropna()
            
            # Need at least 3 variances for Cochran
            if samp_var.size < 3:
                break
                
            r_per_sample = g_work.groupby("Sample")["Replicate"].nunique().mode().iat[0]
            if r_per_sample < 2:
                break
                
            c_res = _cochrans_test(samp_var.values, alpha, df=r_per_sample - 1)
            
            if not c_res.flag:
                break
                
            # Outlier found
            samp_out = samp_var.idxmax()
            variance = samp_var.max()
            
            sample_data = g_work[g_work["Sample"] == samp_out]
            
            outliers.append({
                "Subject": subj,
                "Sample": samp_out,
                "variance": variance,
                "G": c_res.G,
                "G_crit": c_res.G_crit,
                "sample_data": sample_data[["Replicate", "Result"]].to_dict("records")
            })
            
            # Remove from local work set to find next outlier
            g_work = g_work[g_work["Sample"] != samp_out]
            
    return outliers


def _detect_wp_bartlett_outliers(df: pd.DataFrame, alpha: float = 0.05) -> list[dict]:
    """
    Detect QI 10 within-subject variance Bartlett outliers WITHOUT removing them.
    Returns a list of dicts with outlier information for user selection.
    
    Unlike QI 8a which tests replicate variances per sample, QI 10 tests
    within-subject variances (variance of sample means within each subject)
    across subjects.
    
    Each dict contains:
        - Subject: subject ID
        - within_subject_variance: the s²_WP for this subject
        - sample_means: list of {Sample, Mean} dicts
        - chi2: Bartlett chi-squared statistic
        - chi2_crit: critical value
        - p: p-value
    """
    outliers = []
    df_work = df.copy()
    
    while True:
        subj_groups = []
        subj_keys = []
        subj_vars = {}
        
        for subj, g in df_work.groupby("Subject"):
            if g["Sample"].nunique() < 2:
                continue
            sample_means = g.groupby("Sample")["Result"].mean()
            var_wp = sample_means.var(ddof=1)
            if var_wp == 0:
                continue
            subj_groups.append(sample_means.values)
            subj_keys.append(subj)
            subj_vars[subj] = {
                "variance": var_wp,
                "sample_means": sample_means.reset_index().to_dict("records")
            }
        
        if len(subj_groups) < 3:
            break
        
        with np.errstate(all="ignore"):
            bart_stat, bart_p = bartlett(*subj_groups)
        
        k = len(subj_groups)
        chi_crit = chi2.ppf(1 - alpha, k - 1)
        
        if bart_p >= alpha:
            # Homogeneous - done
            break
        
        # Heterogeneous - find worst offender (highest within-subject variance)
        worst_subj = max(subj_vars.keys(), key=lambda s: subj_vars[s]["variance"])
        worst_info = subj_vars[worst_subj]
        
        outliers.append({
            "Subject": worst_subj,
            "within_subject_variance": worst_info["variance"],
            "sample_means": worst_info["sample_means"],
            "chi2": float(bart_stat),
            "chi2_crit": float(chi_crit),
            "p": float(bart_p),
        })
        
        # Remove this subject from working copy to find next outlier
        df_work = df_work[df_work["Subject"] != worst_subj]
    
    return outliers


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


def _detect_reed_outliers(df: pd.DataFrame) -> list[dict]:
    """
    Detect Reed between-subject outlier candidates (QI 8c).
    Returns list of dicts with outlier info (Subject, Mean, Reason).
    """
    outliers = []
    df_work = df.copy()
    
    while True:
        subj_mean = df_work.groupby("Subject")["Result"].mean()
        if subj_mean.size < 3:
            break
            
        ridx = _reed_outlier(subj_mean.values)
        if ridx is None:
            break
            
        culprit = subj_mean.index[ridx]
        mean_val = subj_mean.iloc[ridx]
        
        # Re-calc details for UI context
        values = subj_mean.values
        s_idx = np.argsort(values)
        s_val = values[s_idx]
        rng = s_val[-1] - s_val[0]
        threshold = rng / 3.0
        
        diff_high = s_val[-1] - s_val[-2]
        diff_low = s_val[1] - s_val[0]
        
        reason = "Outlier"
        if (diff_high > threshold) and (subj_mean.index[s_idx[-1]] == culprit):
            reason = f"High-end gap ({diff_high:.2f}) > Range/3 ({threshold:.2f})"
        elif (diff_low > threshold) and (subj_mean.index[s_idx[0]] == culprit):
            reason = f"Low-end gap ({diff_low:.2f}) > Range/3 ({threshold:.2f})"
            
        outliers.append({
            "Subject": culprit,
            "mean": mean_val,
            "range": rng,
            "threshold": threshold,
            "reason": reason
        })
        
        df_work = df_work[df_work["Subject"] != culprit]
        
    return outliers


def _detect_drift_outliers(df: pd.DataFrame) -> list[dict]:
    """
    Detect Steady-State Drift (QI 7) outliers per subject.
    Returns list of dicts with drift info (Subject, Slope, CI, p-value).
    """
    outliers = []
    tinv = lambda p, d: abs(t.ppf(p/2, d))
    
    for subj, g in df.groupby("Subject"):
        if g["Sample"].nunique() <= 2:
            continue
            
        res = linregress(g["Sample"], g["Result"])
        df_denom = len(g) - 2
        
        if df_denom <= 0: continue
        
        if res.pvalue < 0.05:
            ts = tinv(0.05, df_denom)
            ci_low = res.slope - ts * res.stderr
            ci_high = res.slope + ts * res.stderr
            
            outliers.append({
                "Subject": subj,
                "slope": res.slope,
                "p": res.pvalue,
                "ci": (ci_low, ci_high),
                "n_points": len(g)
            })
            
    return outliers



# === Population-wide drift (fixed-effects within estimator) ===
def estimate_population_drift(df: pd.DataFrame,
                              alpha: float = 0.05,
                              collapse_replicates: bool = True) -> dict:
    """
    Estimate a *common within-subject* slope of Result vs Sample across the cohort.
    Approach: subject fixed-effects (within estimator).
      1) Optionally collapse replicates to Subject×Sample means (reduces CVA noise).
      2) Demean both y (Result) and x (Sample) within each subject.
      3) OLS on the demeaned variables → single pooled slope.
      4) Exact t-based CI using df_resid = N - G - 1 (G = #subjects used).

    Returns dict with keys:
      slope, se, t, df, p, ci_low, ci_high, n_subjects, n_points
    """
    if collapse_replicates:
        d = (df.groupby(["Subject", "Sample"], as_index=False)["Result"]
                .mean())
    else:
        # keep replicates; they just add analytical noise
        d = df[["Subject", "Sample", "Result"]].copy()

    # need at least two time points per subject for a slope contribution
    counts = d.groupby("Subject")["Sample"].nunique()
    keep_ids = counts[counts >= 2].index
    d = d[d["Subject"].isin(keep_ids)].copy()

    if d.empty or d["Subject"].nunique() < 2:
        raise ValueError("Insufficient data for population drift (need ≥2 subjects with ≥2 samples).")

    # demean within each subject (within estimator / fixed effects)
    d["y"] = d["Result"] - d.groupby("Subject")["Result"].transform("mean")
    d["x"] = d["Sample"] - d.groupby("Subject")["Sample"].transform("mean")

    # drop rows with zero within-subject x (can happen if all Sample equal)
    d = d[d["x"].abs() > 0].copy()
    if d.empty:
        raise ValueError("No within-subject time variation after demeaning.")

    Sxy = float((d["x"] * d["y"]).sum())
    Sxx = float((d["x"] ** 2).sum())
    slope = Sxy / Sxx

    # residuals & standard error
    resid = d["y"] - slope * d["x"]
    n = len(d)
    g = d["Subject"].nunique()
    df_resid = int(n - g - 1)          # FE: subtract g subject effects + slope
    if df_resid <= 0:
        raise ValueError("Not enough degrees of freedom for pooled slope.")
    sigma2 = float((resid ** 2).sum() / df_resid)
    se = float(np.sqrt(sigma2 / Sxx))

    # t, p, CI
    t_stat = float(slope / se)
    p = float(2 * (1 - t.cdf(abs(t_stat), df_resid)))
    tcrit = float(t.ppf(1 - alpha/2, df_resid))
    ci_low  = float(slope - tcrit * se)
    ci_high = float(slope + tcrit * se)

    return dict(slope=slope, se=se, t=t_stat, df=df_resid, p=p,
                ci_low=ci_low, ci_high=ci_high,
                n_subjects=int(g), n_points=int(n))



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
        flags: dict[str, bool] | None = None,
        rep_outlier_selections: dict | None = None,  # user selections for replicate outliers  (Cochran)
        rep_bartlett_selections: dict | None = None, # user selections for Bartlett outliers
        samp_cochran_selections: dict | None = None, # user selections for Sample level Cochran (QI 8b)
        wp_bartlett_selections: dict | None = None,  # user selections for QI 10 within-subject Bartlett
        reed_selections: dict | None = None,         # user selections for Reed outliers (QI 8c)
        drift_selections: dict | None = None,        # user selections for Steady-state drift (QI 7)
        mapping_warnings: list[str] | None = None,   # warnings from data mapping (e.g., dropped rows)
) -> tuple[pd.DataFrame, list[str]]:
    """
    Apply all pre-analysis quality-improvement (QI) checks and statistical
    filters required by Braga-Panteghini.  Returns a cleaned DataFrame **plus**
    a detailed chronological log (list of strings).
    """
    log: list[str] = []

    # Add mapping warnings at the start of the log (e.g., dropped rows with missing values)
    if mapping_warnings:
        for warn in mapping_warnings:
            log.append(f"Data mapping: {warn}")
    df = df_in.copy()
    transformed = False

    flags = flags or {}                    # safety – empty dict if None
    ON     = lambda name: flags.get(name, True)   # helper: is this step enabled?

    # ════════════════════════════════════════════════════════════════════════
    # 0.  Replicate‑set consistency  (QI 8a – Cochran variance test)
    #     • For every Subject‑×‑Sample pair we calculate the variance of its
    #       replicate measurements.
    #     • Cochran's C identifies the pair whose variance is disproportionately
    #       large.  User selections determine how to handle each outlier.
    # ════════════════════════════════════════════════════════════════════════

    # ── QI 8a – replicate–set Cochran ──────────────────────────────────────────
    if ON("rep_cochran"):
        # If user has provided selections for detected outliers, apply them
        rep_outlier_selections = rep_outlier_selections or {}
        
        # First apply user selections for detected outliers
        for key, action in rep_outlier_selections.items():
            subj, samp = key  # key is (Subject, Sample) tuple
            mask = (df["Subject"] == subj) & (df["Sample"] == samp)
            
            if not mask.any():
                continue  # Already removed or not in data
            
            vals = df.loc[mask, "Result"]
            variance = df.loc[mask, "Result"].var(ddof=1)
            
            if action == "remove_all":
                log.append(
                    f"QI 8a – replicate Cochran outlier removed (user selection) → "
                    f"Subject {subj}, Sample {samp}: "
                    + ", ".join(f"{v:.2f}" for v in vals) +
                    f"  (s² = {variance:.4f})"
                )
                df = df[~mask]
            elif action == "ignore":
                log.append(
                    f"QI 8a – replicate Cochran outlier ignored (user selection) → "
                    f"Subject {subj}, Sample {samp}: "
                    + ", ".join(f"{v:.2f}" for v in vals) +
                    f"  (s² = {variance:.4f}) – kept as-is"
                )
            elif action.startswith("keep_rep_"):
                # Keep only the specified replicate
                rep_to_keep = action.replace("keep_rep_", "")
                # Try to convert to the same type as Replicate column
                try:
                    rep_to_keep = int(float(rep_to_keep))
                except ValueError:
                    pass  # Keep as string if conversion fails
                
                keep_mask = mask & (df["Replicate"] == rep_to_keep)
                remove_mask = mask & (df["Replicate"] != rep_to_keep)
                
                removed_vals = df.loc[remove_mask, "Result"]
                kept_val = df.loc[keep_mask, "Result"]
                
                log.append(
                    f"QI 8a – replicate Cochran outlier: kept Replicate {rep_to_keep} (user selection) → "
                    f"Subject {subj}, Sample {samp}: kept {kept_val.values[0]:.2f}, "
                    f"removed {', '.join(f'{v:.2f}' for v in removed_vals)}"
                )
                df = df[~remove_mask]
        
        # Now run the standard Cochran test on remaining data for any new outliers
        while True:
            rep_var = (df.groupby(["Subject", "Sample"])["Result"]
                        .var(ddof=1).dropna())

            if rep_var.size < 3:
                log.append("QI 8a – Cochran (replicate sets) skipped: <3 replicate-set variances available.")
                break

            R = (df.groupby(["Subject", "Sample"])["Replicate"]
                .nunique().mode().iat[0])

            if R < 2:
                log.append("QI 8a – Cochran (replicate sets) skipped: only 1 replicate/sample.")
                break

            c_res = _cochrans_test(rep_var.values, alpha, df=R - 1)

            log.append(
                f"QI 8a – Cochran (replicate sets): "
                f"G = {c_res.G:.3f}, Gcrit = {c_res.G_crit:.3f}"
                + (" → OUTLIER" if c_res.flag else " → no outlier detected")
            )

            if not c_res.flag:
                break

            # If there are still outliers and no user selection, auto-remove
            subj, samp = rep_var.idxmax()
            if (subj, samp) in rep_outlier_selections:
                # Already handled above, skip
                break
            
            vals = df.loc[(df["Subject"] == subj) & (df["Sample"] == samp), "Result"]
            log.append(
                f"QI 8a – replicate Cochran outlier removed → "
                f"Subject {subj}, Sample {samp}: "
                + ", ".join(f"{v:.2f}" for v in vals) +
                f"  (s² = {rep_var.max():.4f})"
            )
            df = df[~((df["Subject"] == subj) & (df["Sample"] == samp))]
    else:
        log.append("Replicate Cochran (QI 8a) skipped (switch off).")


    # ─────────────────────────────────────────────────────────────────────────
    # ❋ NEW: Bartlett's homogeneity test – analytic variance (replicates)
    # Purpose: Are the replicate variances per Subject×Sample similar? 
    # (Is the analytic imprecison homogeneous across all samples?)
    # ─────────────────────────────────────────────────────────────────────────
    if ON("rep_bartlett"):

        while True:
            valid_groups = []
            group_keys    = []                        # Keep (Subject, Sample) keys
            for (subj, samp), g in df.groupby(["Subject", "Sample"]):
                if g["Result"].size < 2:              # Bartlett needs >=2 measurements
                    continue
                if g["Result"].var(ddof=1) == 0:      # Exclude pairs with s² = 0
                    continue
                valid_groups.append(g["Result"].values)
                group_keys.append((subj, samp))

            if len(valid_groups) >= 3:                # Bartlett needs >=3 groups
                with np.errstate(all="ignore"):       # suppress constant input warnings
                    bart_stat, bart_p = bartlett(*valid_groups)
                rep_msg = "heterogeneous" if bart_p < alpha else "homogeneous"
                k = len(valid_groups)
                chi_crit = chi2.ppf(1 - alpha, k - 1)

                log.append(
                    f"QI 8a – Bartlett test (replicate sets): "
                    f"χ² = {bart_stat:.2f}, χ²crit = {chi_crit:.2f}, "
                    f"p = {bart_p:.4f} → variances {rep_msg}."
                )

                if bart_p < alpha:
                    # Identify outlier variances
                    rep_var = (df.groupby(["Subject", "Sample"])["Result"]
                                .var(ddof=1)
                                .loc[group_keys])     # only those included in test
                    
                    if ON("rep_bartlett_exclusion"):
                        # Iterative exclusion: identify worst offender
                        (bad_subj, bad_samp), max_v = rep_var.idxmax(), rep_var.max()

                        # Check if user has a stored decision for this specific outlier
                        rep_bartlett_selections = rep_bartlett_selections or {}
                        action = rep_bartlett_selections.get((bad_subj, bad_samp), "remove_all")

                        if action == "ignore":
                            log.append(
                                f"  → Bartlett outlier ignored (user selection): "
                                f"Subject {bad_subj}, Sample {bad_samp} (s² = {max_v:.4f}). "
                                "Heterogeneity accepted."
                            )
                            # Stop the loop -> accept current state (heterogeneous)
                            break
                        
                        elif action == "remove_all":
                            log.append(
                                f"  → Removing highest analytical variance outlier: "
                                f"Subject {bad_subj}, Sample {bad_samp} (s² = {max_v:.4f})"
                            )
                            df = df[~((df["Subject"] == bad_subj) & (df["Sample"] == bad_samp))]
                            continue # Loop again

                        elif action.startswith("keep_rep_"):
                            rep_to_keep = action.replace("keep_rep_", "")
                            try:
                                rep_to_keep = int(float(rep_to_keep))
                            except ValueError:
                                pass
                            
                            mask = (df["Subject"] == bad_subj) & (df["Sample"] == bad_samp)
                            keep_mask = mask & (df["Replicate"] == rep_to_keep)
                            kept_val = df.loc[keep_mask, "Result"]
                            
                            log.append(
                                f"  → Bartlett outlier: kept Replicate {rep_to_keep} (user selection) for "
                                f"Subject {bad_subj}, Sample {bad_samp}."
                            )
                            # Keep only that replicate for this sample (effectively variance -> 0 or NaN if <2 reps)
                            # If only 1 replicate remains, it will be excluded from next Bartlett check (size < 2 check)
                            df = df[~mask | keep_mask]
                            continue
                        
                        else:
                            # Fallback default
                            log.append(f"  → Removing outlier (default): Subject {bad_subj}, Sample {bad_samp}")
                            df = df[~((df["Subject"] == bad_subj) & (df["Sample"] == bad_samp))]
                            continue

                    else:
                        # Report only (Default behavior)
                        top_bad = rep_var.sort_values(ascending=False).head(5)
                        for (subj, samp), v in top_bad.items():
                            log.append(
                                f"  High analytical variance → "
                                f"Subject {subj}, Sample {samp}  (s² = {v:.4f})"
                            )
                        break 
                else:
                    # Homogeneous - done
                    break
            else:
                log.append(
                    "QI 8a – Bartlett test (replicate sets) skipped (or stopped): "
                    "fewer than 3 valid replicate groups."
                )
                break
    else:
        log.append("Replicate Bartlett (QI 8a) skipped (switch off).")

    # ── QI 8b – sample‑level Cochran inside each subject ──────────────────
    if ON("samp_cochran"):
        # 1. Apply user selections (if any)
        if samp_cochran_selections:
            for key, action in samp_cochran_selections.items():
                subj, samp = key
                mask = (df["Subject"] == subj) & (df["Sample"] == samp)
                if not mask.any(): continue
                
                vals = df.loc[mask, "Result"]
                variance = df.loc[mask, "Result"].var(ddof=1)
                
                if action == "remove":
                    log.append(f"QI 8b – sample Cochran outlier removed (user selection) → Subject {subj}, Sample {samp} (s²={variance:.4f})")
                    df = df[~mask]
                elif action == "ignore":
                     log.append(f"QI 8b – sample Cochran outlier ignored (user selection) → Subject {subj}, Sample {samp} (s²={variance:.4f})")

        # 2. Auto Mode: Detect and remove all outliers
        if flags.get("outlier_mode") != "manual":
            auto_outliers = _detect_sample_cochran_outliers(df, alpha)
            for out in auto_outliers:
                subj, samp = out["Subject"], out["Sample"]
                vals_str = ", ".join(f"{r['Result']:.2f}" for r in out["sample_data"])
                log.append(f"QI 8b – sample Cochran outlier removed (auto) → Subject {subj}, Sample {samp}: {vals_str} (s²={out['variance']:.4f})")
                df = df[~((df["Subject"] == subj) & (df["Sample"] == samp))]

        # 3. Final Verification / Logging Pass (Check status of retained data)
        for subj, g in df.groupby("Subject"):
            samp_var = g.groupby("Sample")["Result"].var(ddof=1).dropna()
            
            if samp_var.size < 3:
                log.append(f"QI 8b – Cochran skipped (Subject {subj}): <3 variances.")
                continue
                
            r_per_sample = g.groupby("Sample")["Replicate"].nunique().mode().iat[0]
            if r_per_sample < 2:
                log.append(f"QI 8b – Cochran skipped (Subject {subj}): only 1 replicate.")
                continue
                
            c_res = _cochrans_test(samp_var.values, alpha, df=r_per_sample - 1)
            
            status = " → OUTLIER (Retained/Ignored)" if c_res.flag else " → no outlier detected"
            log.append(f"QI 8b – Cochran (Subject {subj}): G = {c_res.G:.3f}, Gcrit = {c_res.G_crit:.3f}{status}")
    else:
        log.append("Sample Cochran (QI 8b) skipped (switch off).")

    # ════════════════════════════════════════════════════════════════════════
    # QI 10 – Within-subject Bartlett (variance homogeneity across subjects)
    # ════════════════════════════════════════════════════════════════════════
    if ON("wp_bartlett"):
        wp_bartlett_selections = wp_bartlett_selections or {}
        
        while True:
            subj_groups = []
            subj_keys = []
            subj_vars = {}
            
            for subj, g in df.groupby("Subject"):
                if g["Sample"].nunique() < 2:
                    continue
                sample_means = g.groupby("Sample")["Result"].mean()
                var_wp = sample_means.var(ddof=1)
                if var_wp == 0:
                    continue
                subj_groups.append(sample_means.values)
                subj_keys.append(subj)
                subj_vars[subj] = var_wp
            
            if len(subj_groups) < 3:
                log.append(
                    "QI 10 – Bartlett test (within‑subject variances) skipped: "
                    "fewer than 3 eligible subjects."
                )
                break
            
            with np.errstate(all="ignore"):
                bart_stat, bart_p = bartlett(*subj_groups)
            
            k = len(subj_groups)
            chi_crit = chi2.ppf(1 - alpha, k - 1)
            wp_msg = "heterogeneous" if bart_p < alpha else "homogeneous"
            
            log.append(
                f"QI 10 – Bartlett test (within‑subject variances): "
                f"χ² = {bart_stat:.2f}, χ²crit = {chi_crit:.2f}, "
                f"p = {bart_p:.4f} → variances {wp_msg}."
            )
            
            if bart_p >= alpha:
                # Homogeneous - done
                break
            
            # Heterogeneous - check if exclusion is enabled
            if ON("wp_bartlett_exclusion"):
                # Find worst offender
                worst_subj = max(subj_vars.keys(), key=lambda s: subj_vars[s])
                max_var = subj_vars[worst_subj]
                
                # Check user selection
                action = wp_bartlett_selections.get(worst_subj, "remove")
                
                if action == "ignore":
                    log.append(
                        f"  → QI 10 outlier ignored (user selection): "
                        f"Subject {worst_subj} (s²_WP = {max_var:.4f}). "
                        "Heterogeneity accepted."
                    )
                    break
                elif action == "remove":
                    log.append(
                        f"  → Removing highest within-subject variance: "
                        f"Subject {worst_subj} (s²_WP = {max_var:.4f})"
                    )
                    df = df[df["Subject"] != worst_subj]
                    continue
                else:
                    # Fallback - remove
                    log.append(f"  → Removing outlier (default): Subject {worst_subj}")
                    df = df[df["Subject"] != worst_subj]
                    continue
            else:
                # Report only (no exclusion enabled)
                for subj, v in sorted(subj_vars.items(), key=lambda x: x[1], reverse=True)[:5]:
                    log.append(
                        f"  High within‑subject variance → Subject {subj}  (s² = {v:.4f})"
                    )
                break
    else:
        log.append("Within‑subject Bartlett (QI 10) skipped (switch off).")

    # ── QI 8c – between‑subject Reed test only ────────────────────────────────
    if ON("reed"):
        # 1. Apply user selections (if any)
        if reed_selections:
            for subj, action in reed_selections.items():
                if action == "remove":
                    log.append(f"QI 8c – Reed outlier removed (user selection) → Subject {subj}")
                    df = df[df["Subject"] != subj]
                elif action == "ignore":
                    log.append(f"QI 8c – Reed outlier ignored (user selection) → Subject {subj}")

        # 2. Auto Mode
        if flags.get("outlier_mode") != "manual":
            auto_reed = _detect_reed_outliers(df)
            for out in auto_reed:
                subj = out["Subject"]
                log.append(f"QI 8c – Reed mean outlier removed (auto) → Subject {subj} (mean={out['mean']:.4f})")
                df = df[df["Subject"] != subj]
        
        # 3. Final Verification
        subj_mean = df.groupby("Subject")["Result"].mean()
        if subj_mean.size >= 3:
            ridx = _reed_outlier(subj_mean.values)
            if ridx is None:
                log.append("QI 8c – Reed test passed: no between‑subject outliers.")
            else:
                culprit = subj_mean.index[ridx]
                log.append(f"QI 8c – Reed outlier identified (retained/ignored): Subject {culprit}.")
    else:
        log.append("Reed between‑subject (QI 8c) skipped (switch off).")


    # ════════════════════════════════════════════════════════════════════════
    # 0 b.  Steady‑state trend test (QI 7) – slope ± 95 % CI
    # ════════════════════════════════════════════════════════════════════════
    if ON("drift"):
        # 1. Apply user selections (if any)
        if drift_selections:
            for subj, action in drift_selections.items():
                if action == "remove":
                    log.append(f"QI 7 – Drift outlier removed (user selection) → Subject {subj}")
                    df = df[df["Subject"] != subj]
                elif action == "ignore":
                    log.append(f"QI 7 – Drift outlier ignored (user selection) → Subject {subj}")

        # 2. Auto Mode
        if flags.get("outlier_mode") != "manual":
            auto_drift = _detect_drift_outliers(df)
            for out in auto_drift:
                subj = out["Subject"]
                log.append(
                    f"QI 7 – temporal drift removed (auto): Subject {subj} "
                    f"slope={out['slope']:.3g} (p={out['p']:.3g})"
                )
                df = df[df["Subject"] != subj]

        # 3. Final Verification
        final_drifts = _detect_drift_outliers(df)
        if not final_drifts:
            log.append("QI 7 – Steady-state drift test passed: no trend outliers.")
        else:
            log.append(f"QI 7 – Drift outliers identified (retained/ignored): {len(final_drifts)} subjects.")
            
    else:
        log.append("Steady‑state drift (QI 7) skipped (switch off).")


    # ─────────────────────────────────────────────────────────────
    # QI 7 (population) – pooled within-subject drift across cohort
    # ─────────────────────────────────────────────────────────────
    if ON("pop_drift"):
        try:
            pop = estimate_population_drift(df, alpha=alpha, collapse_replicates=True)
            # report slope in “Result unit per one Sample step”
            unit = st.session_state.get("result_unit", "").strip() if "result_unit" in st.session_state else ""
            unit_txt = f" {unit}" if unit else ""
            log.append(
                "QI 7 – population drift: "
                f"slope = {pop['slope']:+.3g}{unit_txt}/sample "
                f"(95% CI {pop['ci_low']:.3g} to {pop['ci_high']:.3g}, "
                f"t = {pop['t']:.2f}, df = {pop['df']}, p = {pop['p']:.3g}; "
                f"{pop['n_subjects']} subjects; {pop['n_points']} points). "
                "No exclusions are performed at this step."
            )

        except Exception as e:
            log.append(f"QI 7 – population drift estimation skipped: {e}")
    else:
        log.append("Population drift (QI 7, pooled) skipped (switch off).")


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
                p_sw_log = shapiro(means)[1]
                if p_sw_log <= normal_p:
                    log.append(f"WARNING: Normality could not be achieved even after log-transform (p={p_sw_log:.3g}) – (Proceeding with log-data, but interpretation requires caution).")
            elif p_sw <= normal_p:
                log.append(f"WARNING: Normality could not be achieved (p={p_sw:.3g}) – (Proceeding, but interpretation requires caution).")
    else:
        log.append("Normality checks / log‑transform (QI 9) skipped (switch off).")  



    # ════════════════════════════════════════════════════════════════════════
    # 4.  Force a perfectly balanced design for the ANOVA
    # ════════════════════════════════════════════════════════════════════════
    # ════════════════════════════════════════════════════════════════════════
    # 4.  Force a perfectly balanced design (equal S & R)
    #     – required for the closed-form two-level ANOVA that follows
    #     NOTE: We run sample-level (replicate) balance FIRST, then subject-level,
    #           so that if a sample is removed due to unequal replicates, the
    #           subject-level check will catch subjects with fewer samples.
    # ════════════════════════════════════════════════════════════════════════
    # ═════════ 4.  Force a perfectly balanced design ══════════════════════
    if enforce_balance:                                      # NEW GUARD

        # 4-a ▸ SAMPLE-LEVEL balance (replicate count) — RUN FIRST ─────────────
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
            bad_idx = df.set_index(["Subject", "Sample"]).index.isin(
                bad_pairs.set_index(["Subject", "Sample"]).index)
            df = df[~bad_idx]

        # 4-b ▸ SUBJECT-LEVEL balance (sample count) — RUN AFTER ───────────────
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

def _cvg_ln_balanced(clean_df: pd.DataFrame, alpha: float
) -> tuple[float, tuple[float, float]]:
    df = clean_df.copy()
    if (df["Result"] <= 0).any():
        raise ValueError("CVG (ln) requires all results > 0.")
    df["Z"] = np.log(df["Result"])

    I = df["Subject"].nunique()
    S = df.groupby("Subject")["Sample"].nunique().iloc[0]
    R = df.groupby(["Subject", "Sample"])["Replicate"].nunique().iloc[0]

    subj_mean = df.groupby("Subject")["Z"].mean()
    samp_mean = df.groupby(["Subject", "Sample"])["Z"].mean()
    grand     = df["Z"].mean()

    ss_bp = S*R * ((subj_mean - grand) ** 2).sum()
    ss_wp = R   * ((samp_mean - samp_mean.index.get_level_values(0).map(subj_mean)) ** 2).sum()

    df_bp = I - 1
    df_wp = I * (S - 1)

    ms_bp = ss_bp / df_bp
    ms_wp = ss_wp / df_wp

    var_bp_ln = max((ms_bp - ms_wp) / (S * R), 0.0)
    cvg       = float(np.sqrt(np.exp(var_bp_ln) - 1.0) * 100)

    ms_lo = (df_bp * ms_bp) / chi2.ppf(1 - alpha/2, df_bp)
    ms_hi = (df_bp * ms_bp) / chi2.ppf(alpha/2,     df_bp)

    var_lo = max((ms_lo - ms_wp) / (S * R), 0.0)
    var_hi = max((ms_hi - ms_wp) / (S * R), 0.0)

    ci = (float(np.sqrt(np.exp(var_lo) - 1.0) * 100),
          float(np.sqrt(np.exp(var_hi) - 1.0) * 100))
    return cvg, ci


def _cvg_ln_unbalanced(clean_df: pd.DataFrame, alpha: float
) -> tuple[float, tuple[float, float]]:
    df = clean_df.copy()
    if (df["Result"] <= 0).any():
        raise ValueError("CVG (ln) requires all results > 0.")
    df["Z"] = np.log(df["Result"])

    I    = df["Subject"].nunique()
    S_i  = df.groupby("Subject")["Sample"].nunique()
    r_ij = df.groupby(["Subject","Sample"])["Replicate"].nunique()
    n_i  = r_ij.groupby("Subject").sum()
    nbar = float(n_i.mean())                 # = S̄·r̄

    subj_mean = df.groupby("Subject")["Z"].mean()
    samp_mean = df.groupby(["Subject","Sample"])["Z"].mean()
    grand     = df["Z"].mean()

    ss_bp = (n_i * (subj_mean - grand)**2).sum()
    ss_wp = ((samp_mean - samp_mean.index.get_level_values(0).map(subj_mean))**2 * r_ij).sum()

    df_bp = I - 1
    df_wp = (S_i - 1).sum()

    ms_bp = ss_bp / df_bp
    ms_wp = ss_wp / df_wp

    var_bp_ln = max((ms_bp - ms_wp) / nbar, 0.0)
    cvg       = float(np.sqrt(np.exp(var_bp_ln) - 1.0) * 100)

    ms_lo = (df_bp * ms_bp) / chi2.ppf(1 - alpha/2, df_bp)
    ms_hi = (df_bp * ms_bp) / chi2.ppf(alpha/2,     df_bp)

    var_lo = max((ms_lo - ms_wp) / nbar, 0.0)
    var_hi = max((ms_hi - ms_wp) / nbar, 0.0)

    ci = (float(np.sqrt(np.exp(var_lo) - 1.0) * 100),
          float(np.sqrt(np.exp(var_hi) - 1.0) * 100))
    return cvg, ci


# === Lognormal RCV helper (Fokkema et al.) ===
Z_BIDIR_95 = 1.64  # one-sided 95% by default

def rcv_lognormal(cv_a_pct: float, cv_i_pct: float, z: float = Z_BIDIR_95) -> tuple[float, float]:
    """
    Return (RCV_down %, RCV_up %) using the lognormal/asymmetric method.

    σ_A² = ln(1 + (CVa/100)²),  σ_I² = ln(1 + (CVi/100)²)
    SD   = sqrt( 2 * (σ_A² + σ_I²) )
    RCV↑ = (exp(+z * SD) - 1) * 100
    RCV↓ = (1 - exp(-z * SD)) * 100
    """
    cv_a = cv_a_pct / 100.0
    cv_i = cv_i_pct / 100.0
    sigma2 = np.log1p(cv_a**2) + np.log1p(cv_i**2)  # ln(1 + cv^2)
    SD = np.sqrt(2.0 * sigma2)
    rcv_up = (np.exp(z * SD) - 1.0) * 100.0
    rcv_down = (1.0 - np.exp(-z * SD)) * 100.0
    return rcv_down, rcv_up

# === APS from BV (imprecision, bias, MAU) ===
def build_aps_table(cvi_pct: float, cvg_pct: float, k: float = 2.0) -> pd.DataFrame:
    """
    Return APS (Minimum / Desirable / Optimal) as a DataFrame.

    Formulas (all in %):
      CVa_spec  = f × CVI
      Bias_spec = g × sqrt(CVI² + CVG²)
      MAU_spec  = k × f × CVI        # expanded allowable measurement uncertainty

    Factors by level:
      Optimal:   f=0.25, g=0.125
      Desirable: f=0.50, g=0.25
      Minimum:   f=0.75, g=0.375

    Abbreviations: CVa (analytical imprecision), CVI (within-subject BV),
                   CVG (between-subject BV), MAU (expanded allowable measurement uncertainty).
    """
    term = float(np.sqrt(cvi_pct**2 + cvg_pct**2))
    levels = [
        ("Minimum",   0.75,  0.375),
        ("Desirable", 0.50,  0.25),
        ("Optimal",   0.25,  0.125),
    ]
    rows = []
    for name, f_imp, g_bias in levels:
        rows.append({
            "Specification": name,
            "CVa":  f_imp * cvi_pct,
            "Bias": g_bias * term,
            "MAU (k=2)": k * f_imp * cvi_pct,  # k fixed to 2
        })
    return pd.DataFrame(rows)


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

    cvg_ln, ci_cvg_ln = _cvg_ln_unbalanced(df, alpha)
    cv_G     = cvg_ln
    ci_cv_G  = ci_cvg_ln
                               
    rcv_down, rcv_up = rcv_lognormal(cv_A, cv_I, z=Z_BIDIR_95)


    return dict(var_A=var_A, var_WP=var_WP, var_BP=var_BP,
                cv_A=cv_A, cv_I=cv_I, cv_G=cv_G,
                ci_cv_A=ci_cv_A, ci_cv_I=ci_cv_I, ci_cv_G=ci_cv_G,
                rcv_down=rcv_down, rcv_up=rcv_up, grand=grand,
                df_bp=df_bp, df_wp=df_wp, df_a=df_a)


# ——————————————————————————————————————————————————————————————
# 4.  Core calculation routine (closed‑form Røraas algebra)
# ——————————————————————————————————————————————————————————————

def calculate_bv(df: pd.DataFrame, alpha: float = 0.05,
                 use_cv_anova: bool = False,
                 enforce_balance: bool = True,
                 rep_outlier_selections: dict | None = None,  # user selections for replicate outliers (Cochran)
                 rep_bartlett_selections: dict | None = None, # user selections for Bartlett outliers
                 samp_cochran_selections: dict | None = None, # user selections for Sample level Cochran (QI 8b)
                 wp_bartlett_selections: dict | None = None,  # user selections for QI 10 within-subject Bartlett
                 reed_selections: dict | None = None,         # user selections for Reed outliers (QI 8c)
                 drift_selections: dict | None = None,        # user selections for Steady-state drift (QI 7)
                 mapping_warnings: list[str] | None = None,   # warnings from data mapping (e.g., dropped rows)
) -> BVResult:
    """Compute balanced‑design variance components & CVs.

    Parameters
    ----------
    df : DataFrame
        Must contain at least the four columns: Subject, Sample, Replicate, Result.
        Extra columns are ignored. Capitalisation is normalised.
    alpha : float
        Significance level for the CI on CV_I (0.05 → 95 % CI).
    rep_outlier_selections : dict | None
        User selections for how to handle detected replicate outliers (Cochran).
        Keys are (Subject, Sample) tuples, values are actions like "remove_all", "ignore", or "keep_rep_N".
    rep_bartlett_selections : dict | None
        User selections for how to handle Bartlett outliers.
        Keys are (Subject, Sample) tuples, values are actions like "remove_all", "ignore", or "keep_rep_N".
    samp_cochran_selections : dict | None
        User selections for Sample level Cochran outliers.
        Keys are (Subject, Sample) tuples, values are "remove" or "ignore".
    wp_bartlett_selections : dict | None
        User selections for QI 10 within-subject Bartlett outliers.
        Keys are Subject IDs, values are "remove" or "ignore".
    reed_selections : dict | None
        User selections for Reed outliers (QI 8c).
        Keys are Subject IDs, values are "remove" or "ignore".
    drift_selections : dict | None
        User selections for Steady-state drift (QI 7).
        Keys are Subject IDs, values are "remove" or "ignore".
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
            df,
            alpha=alpha,
            enforce_balance=enforce_balance,
            flags=st.session_state.get("preproc_flags"),
            rep_outlier_selections=rep_outlier_selections,
            rep_bartlett_selections=rep_bartlett_selections,
            samp_cochran_selections=samp_cochran_selections,
            wp_bartlett_selections=wp_bartlett_selections,
            reed_selections=reed_selections,
            drift_selections=drift_selections,
            mapping_warnings=mapping_warnings or [],
        )
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
            rcv_95_down = ub["rcv_down"],
            rcv_95_up   = ub["rcv_up"],
            preprocess_log = pp_log,
            clean_df=df.copy(),   # NEW
            cv_I_cv_anova = cv_anova,               
            ci_cv_I_cv_anova = ci_cv_anova,   
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
        # --- CVG: ln-ANOVA + back-transform ---
        cv_G, ci_cv_G = _cvg_ln_balanced(df, alpha)
        pp_log.append("CVG estimated on ln scale and back-transformed.")

        # 3.11 reference‑change value (two‑sided, 95 %)
        rcv_down, rcv_up = rcv_lognormal(cv_A, cv_I, z=Z_BIDIR_95)
        
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
            ci_cv_I, ci_cv_G, rcv_down, rcv_up,
            preprocess_log = pp_log,
            clean_df=df.copy(),   # NEW
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
    # read unit from session_state if present
    unit = st.session_state.get("result_unit", "").strip() if "result_unit" in st.session_state else ""
    measurand = st.session_state.get("measurand_name", "").strip() if "measurand_name" in st.session_state else "Result"
    if not measurand: measurand = "Result"
    y_title = f"{measurand}" + (f" ({unit})" if unit else "")

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
            title=y_title,
            zeroline=False,
            showgrid=False,
            ticks="outside",
            ticklen=4,
            tickcolor="#bfc1c4"
        ),
        font=dict(size=11),
        title=dict(
            text="Per-subject mean (range)",
            x=0.5, xanchor="center",
            font=dict(size=14, color="#2a3f5f")
        ),
        hoverlabel=dict(bgcolor="white", font_size=11),
        showlegend=False
    )

    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Gender-stratified subject range plot
# ─────────────────────────────────────────────────────────────────────────────
def plot_gender_stratified_ranges(clean_df: pd.DataFrame, gender_map: dict[str, str]) -> go.Figure:
    """
    Plot per-subject ranges (min-max) and means, stratified by Gender (Female/Male).
    Includes:
      - Group Mean (solid line)
      - 95% CI of the Mean (shaded band)
      - 5th and 95th Percentiles (dashed lines)
    """
    # 1. Merge gender info
    df = clean_df.copy()
    df["Gender"] = df["Subject"].map(gender_map)
    
    # Filter to known genders and required columns
    df = df[df["Gender"].isin(["Male", "Female"])].copy()
    if df.empty:
        return go.Figure().update_layout(title="No gender-mapped data available for plot")

    # 2. Aggregations per subject
    subj_agg = (df.groupby(["Gender", "Subject"])["Result"]
                  .agg(mean="mean", min="min", max="max")
                  .reset_index())
    
    # Sort: Female first, then Male. Within gender: sort by Subject ID (or mean?)
    subj_agg = subj_agg.sort_values(by=["Gender", "Subject"], ascending=[True, True])
    
    # 3. Calculate Group Statistics
    stats = {}
    for g in ["Female", "Male"]:
        g_data = subj_agg[subj_agg["Gender"] == g]
        if g_data.empty:
            continue
            
        group_mean = g_data["mean"].mean()
        # 95% CI of the Mean: Mean +/- t * SE
        n = len(g_data)
        se = g_data["mean"].std(ddof=1) / np.sqrt(n) if n > 1 else 0
        t_crit = t.ppf(0.975, n-1) if n > 1 else 0
        ci_lower = group_mean - t_crit * se
        ci_upper = group_mean + t_crit * se
        
        # Percentiles (5th and 95th) of the population (subject means)
        p05 = g_data["mean"].quantile(0.05)
        p95 = g_data["mean"].quantile(0.95)
        
        stats[g] = {
            "mean": group_mean,
            "ci_low": ci_lower,
            "ci_high": ci_upper,
            "p05": p05,
            "p95": p95,
            "count": n
        }

    # 4. Plotting
    fig = go.Figure()
    
    subj_agg["x_label"] = subj_agg["Subject"].astype(str).apply(lambda x: f"Subject {x}")
    
    current_x = 0
    boundary_x = None
    
    configs = {
        "Female": {"color": "#d62728", "name": "Female"}, # Red
        "Male":   {"color": "#1f77b4", "name": "Male"}    # Blue
    }
    
    for gender in ["Female", "Male"]:
        g_df = subj_agg[subj_agg["Gender"] == gender]
        if g_df.empty:
            continue
            
        # Add traces for subjects
        fig.add_trace(go.Scatter(
            x=g_df["x_label"], 
            y=g_df["mean"],
            mode='markers',
            name=f"{gender} Subjects",
            marker=dict(color=configs[gender]["color"], size=6),
            error_y=dict(
                type='data',
                symmetric=False,
                array=g_df["max"] - g_df["mean"],
                arrayminus=g_df["mean"] - g_df["min"],
                color=configs[gender]["color"],
                thickness=1,
                width=0
            ),
            showlegend=False,
            hovertemplate="Subject: %{x}<br>Mean: %{y:.2f}<extra></extra>"
        ))
        
        # Add Reference Lines
        s = stats.get(gender)
        if s:
            x0 = current_x - 0.4
            x1 = current_x + s["count"] - 0.6
            
            # Mean (Solid)
            fig.add_shape(type="line", x0=x0, x1=x1, y0=s["mean"], y1=s["mean"],
                          line=dict(color=configs[gender]["color"], width=2, dash="solid"))
            
            # CI Band (Shaded Rectangle)
            fig.add_shape(type="rect", x0=x0, x1=x1, y0=s["ci_low"], y1=s["ci_high"],
                          fillcolor=configs[gender]["color"], opacity=0.15, line_width=0)
            
            # Percentiles (Dashed)
            fig.add_shape(type="line", x0=x0, x1=x1, y0=s["p05"], y1=s["p05"],
                          line=dict(color=configs[gender]["color"], width=1, dash="dash"))
            fig.add_shape(type="line", x0=x0, x1=x1, y0=s["p95"], y1=s["p95"],
                          line=dict(color=configs[gender]["color"], width=1, dash="dash"))
        
        current_x += len(g_df)
        
        if gender == "Female" and "Male" in stats:
            boundary_x = current_x - 0.5

    # Vertical Divider Line
    if boundary_x is not None:
        fig.add_vline(x=boundary_x, line_width=1, line_dash="dash", line_color="black")

    # Layout
    unit = st.session_state.get("result_unit", "").strip() if "result_unit" in st.session_state else ""
    measurand = st.session_state.get("measurand_name", "").strip() if "measurand_name" in st.session_state else "Result"
    if not measurand: measurand = "Result"
    y_title = f"{measurand}" + (f" ({unit})" if unit else "")

    fig.update_layout(
        title="Gender-stratified Per-subject Distribution",
        template="simple_white",
        yaxis_title=y_title,
        xaxis_title="",
        xaxis=dict(
            tickangle=-45,
            showgrid=True,
            gridcolor="#f0f0f0"
        ),
        yaxis=dict(
            showgrid=True,
            gridcolor="#f0f0f0"
        ),
        showlegend=True,
        height=500
    )
    
    # Manual Legend items
    for g in ["Female", "Male"]:
        if g in stats:
            # Group Mean
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode='lines',
                line=dict(color=configs[g]["color"], dash="solid", width=2),
                name=f"{g} Mean"
            ))
            # 95% CI of Mean (Shaded Region)
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode='markers',
                marker=dict(symbol="square", color=configs[g]["color"], opacity=0.3, size=10),
                name=f"{g} Mean 95% CI"
            ))
            # 5th and 95th Percentiles
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode='lines',
                line=dict(color=configs[g]["color"], dash="dash", width=1),
                name=f"{g} 5th-95th Percentile"
            ))
            
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

# ─────────────────────────────────────────────────────────────────────────────
# Gender-based analysis helpers
# ─────────────────────────────────────────────────────────────────────────────
def _enforce_subject_level_gender(df_gender: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """
    Ensure each Subject has exactly one gender label across all its rows.
    Drop subjects with missing/ambiguous/mixed gender.
    """
    glog = []

    uniq = df_gender.groupby("Subject")["Gender"].apply(lambda s: list(pd.unique(s.dropna())))
    good = uniq[uniq.apply(len) == 1]
    bad  = uniq[uniq.apply(len) != 1]

    if len(bad):
        glog.append(f"Gender mapping: excluded {len(bad)} subject(s) with missing/ambiguous/mixed gender.")

    keep_ids = set(good.index)
    out = df_gender[df_gender["Subject"].isin(keep_ids)].copy()
    subj_to_gender = {sid: genders[0] for sid, genders in good.items()}
    out["Gender"] = out["Subject"].map(subj_to_gender)
    return out, glog

def ci_overlaps(ci1: tuple[float, float], ci2: tuple[float, float]) -> bool:
    """True if two closed intervals overlap (touching counts as overlap)."""
    if ci1 is None or ci2 is None:
        return True
    lo = max(float(ci1[0]), float(ci2[0]))
    hi = min(float(ci1[1]), float(ci2[1]))
    return lo <= hi

def _derive_gender_single_column(raw_df: pd.DataFrame,
                                gender_col: str,
                                male_value,
                                female_value) -> pd.Series:
    """Map a single gender column into canonical labels: 'Male'/'Female'/NaN."""
    s = raw_df[gender_col]
    g = pd.Series(index=raw_df.index, dtype="object")
    g.loc[s == male_value] = "Male"
    g.loc[s == female_value] = "Female"
    return g

def _derive_gender_two_indicators(raw_df: pd.DataFrame,
                                 male_col: str,
                                 female_col: str) -> pd.Series:
    """
    Map two indicator columns into canonical labels.
    Truthiness rule: nonzero/True => indicator set.
    Ambiguous (both True or both False/NA) => NaN.
    """
    m = raw_df[male_col].fillna(0)
    f = raw_df[female_col].fillna(0)

    m_true = m.astype(bool)
    f_true = f.astype(bool)

    g = pd.Series(index=raw_df.index, dtype="object")
    g.loc[m_true & ~f_true] = "Male"
    g.loc[f_true & ~m_true] = "Female"
    # ambiguous stays NaN
    return g

def _pick_cvi_ci(res: BVResult, use_cv_anova_flag: bool) -> tuple[float, tuple[float, float]]:
    """
    Decide which CVI value/CI is “the CVI” for overlap decisions.
    If CV-ANOVA is selected AND available, use that; otherwise use standard ANOVA.
    """
    if use_cv_anova_flag and (res.cv_I_cv_anova is not None) and (res.ci_cv_I_cv_anova is not None):
        return float(res.cv_I_cv_anova), tuple(res.ci_cv_I_cv_anova)
    return float(res.cv_I), tuple(res.ci_cv_I)

def _render_gender_overlap_recommendations(res_m: BVResult, res_f: BVResult, use_cv_anova_flag: bool):
    # CVI overlap check
    _, cvi_ci_m = _pick_cvi_ci(res_m, use_cv_anova_flag)
    _, cvi_ci_f = _pick_cvi_ci(res_f, use_cv_anova_flag)

    if not ci_overlaps(cvi_ci_m, cvi_ci_f):
        st.warning(
            "✅ **Gender-based recommendation (CVI):** CVI 95% CIs **do NOT overlap** between genders. "
            "**Report CVI separately** for Male and Female."
        )
    else:
        st.info(
            "CVI 95% CIs **overlap** between genders → separate CVI reporting is not required."
        )

    # Mean overlap check (drives CVG recommendation per your requirement)
    if not ci_overlaps(res_m.ci_mean, res_f.ci_mean):
        st.warning(
            "✅ **Gender-based recommendation (CVG):** Mean concentration 95% CIs **do NOT overlap** between genders. "
            "**Report CVG separately** for Male and Female."
        )
    else:
        st.info(
            "Mean concentration 95% CIs **overlap** between genders → separate CVG reporting is not required."
        )


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
    st.sidebar.subheader("Outlier Exclusion")


    # Mode selector: Exclude all vs Manual review
    outlier_mode = st.sidebar.radio(
        "Outlier handling mode:",
        ["🔄 Exclude all outliers (automatic)", "🔍 Manually review outliers"],
        index=0,
        help="Choose whether to automatically exclude all detected outliers or review each one individually."
    )
    st.session_state["outlier_mode"] = "auto" if outlier_mode.startswith("🔄") else "manual"
    
    # Outlier exclusion checkboxes (in specified order)
    OUTLIER_OPTS = {
        "Replicate level Cochran (QI 8a)"     : ("rep_cochran",   True),
        "Replicate level Bartlett (QI 8a)"    : ("rep_bartlett",  True),
        "Sample level Cochran (QI 8b)"        : ("samp_cochran",  True),
        "Within-subject Bartlett (QI 10)"     : ("wp_bartlett",   True),
        "Reed between-subject (QI 8c)"        : ("reed",          True),
        "Steady-state drift (QI 7)"           : ("drift",         True),
    }
    
    preproc_flags = {}
    for label, (key, default) in OUTLIER_OPTS.items():
        preproc_flags[key] = st.sidebar.checkbox(label, value=default)
    
    # QI 7: If drift is NOT selected, ask if trend analysis is important
    if not preproc_flags.get("drift", True):
        qi7_importance = st.sidebar.radio(
            "Is trend analysis important for this measurand?",
            options=["Unlikely important (Grade B)", "May be important (Grade C)"],
            index=1,  # Default to Grade C (conservative)
            help="BIVAC QI 7: If trend analysis is not performed, grade depends on importance for the measurand."
        )
        st.session_state["qi7_importance"] = "B" if "Grade B" in qi7_importance else "C"
    else:
        st.session_state["qi7_importance"] = "A"  # Drift test performed = Grade A
    
    # Further Analysis Section
    st.sidebar.subheader("Further Analysis")
    ANALYSIS_OPTS = {
        "Normality checks / log-transform (QI 9)" : ("normality", True),
        "Population drift (QI 7)"                 : ("pop_drift", False),
    }
    for label, (key, default) in ANALYSIS_OPTS.items():
        preproc_flags[key] = st.sidebar.checkbox(label, value=default)
    
    # Implicitly enable exclusion logic if the parent check is ON
    preproc_flags["rep_bartlett_exclusion"] = preproc_flags.get("rep_bartlett", False)
    preproc_flags["wp_bartlett_exclusion"] = preproc_flags.get("wp_bartlett", False)
    
    # Include QI 7 importance setting
    preproc_flags["qi7_importance"] = st.session_state.get("qi7_importance", "C")
    
    preproc_flags["outlier_mode"] = st.session_state.get("outlier_mode", "auto")
    st.session_state["preproc_flags"] = preproc_flags

    # Analysis Options Section
    st.sidebar.subheader("Analysis Options")
    enforce_balance = st.sidebar.checkbox(
        "Enforce balanced crossed design",
        value=True,
    )
    st.session_state["enforce_balance"] = enforce_balance

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


# 6.3 tabs – add an Instructions tab
instr_tab, upload_tab, entry_tab = st.tabs(["**Instructions**", "**Upload**", "**Manual Entry**"])

with instr_tab:
    st.subheader("How this app works (quick guide)")
    st.markdown(
        """
1. **Prepare your data** as a tidy table with four columns:
   - **Subject** – participant ID (e.g., `P01`, `MRN123`).
   - **Sample / time-point** – visit/order **as integers**: `1, 2, 3…`.
   - **Replicate** – within-sample repeat **as integers**: `1, 2…`.
   - **Result / value** – the numeric measurement (one unit).
2. **Upload** a CSV/XLSX or use **Manual Entry**.
3. **Map / confirm columns** in the expander and click **Save mapping**.
4. (Optional) Adjust **Pre-processing steps** and **Analysis options** in the sidebar.
5. Click **Calculate** to get CVₐ, CVᵢ, CVg, 95% CIs, RCV, per-subject plot, and a BIVAC checklist.
        """
    )

    st.divider()
    c1, c2 = st.columns([0.55, 0.45], gap="large")

    with c1:
        st.markdown("### Data format – minimum rules")
        st.markdown(
            """
- **Subject** can be text or numbers (blanks discouraged).
- **Sample** and **Replicate** **must be numeric**; non-numeric entries will error.
- **Result** **must be numeric**; non-numeric rows are rejected during mapping.
- Balanced designs (**same S and R for all subjects**) enable closed-form ANOVA.
- Unbalanced is allowed if you untick **“Enforce balanced crossed design”** (sidebar).
            """
        )
        st.markdown("**Example structure (first rows of the built-in template):**")
        st.dataframe(_template_df.head(8), use_container_width=True, hide_index=True)
       
        st.markdown("### Outputs you get")
        st.markdown(
            """
- **Key metrics**: Mean, **CVₐ**, **CVᵢ**, **CVg**, 95% CIs, and **RCV (95%)**.
- **Per-subject mean ± range** plot.
- **BIVAC checklist (QI 1–14)** with export to Excel.
- **Pre-processing log** with highlighted actions/outliers.
            """
        )

    with c2:
        st.markdown("### What the pre-processing switches do")
        st.markdown(
            """
- **Replicate Cochran (QI 8a)**: flags replicate-set variance outliers.
- **Replicate Bartlett (QI 8a)**: checks homogeneity of replicate variances.
- **Sample Cochran (QI 8b)**: flags high-variance samples within each subject.
- **Reed between-subject (QI 8c)**: removes extreme subject means.
- **Steady-state drift (QI 7)**: removes subjects with significant time trend.
- **Normality checks / log-transform (QI 9)**: Shapiro–Wilk; log if needed.
- **Within-subject Bartlett (QI 10)**: homogeneity of within-subject variances.
- **Enforce balanced crossed design**: keeps only subjects/samples that fit the modal S and R.
- **Estimate CVI with CV-ANOVA**: optional CVᵢ via normalization approach.
            """
        )

        # st.markdown("### Recommended study design")
        # st.info("**≥20 subjects**, **≥3 samples/subject**, **≥2 replicates/sample** is a good starting point.", icon="✅")
    st.info("**When using BIVAC, please cite the following reference:** *Aarsand AK, Røraas T, Fernandez-Calle P, Ricos C, Díaz-Garzón J, Jonker N, Perich C, González-Lao E, Carobene A, Minchinela J, Coşkun A, Simón M, Álvarez V, Bartlett WA, Fernández-Fernández P, Boned B, Braga F, Corte Z, Aslan B, Sandberg S; European Federation of Clinical Chemistry and Laboratory Medicine Working Group on Biological Variation and Task and Finish Group for the Biological Variation Database. The Biological Variation Data Critical Appraisal Checklist: A Standard for Evaluating Studies on Biological Variation. Clin Chem. 2018 Mar;64(3):501-514. doi: 10.1373/clinchem.2017.281808. Epub 2017 Dec 8. PMID: 29222339.*")
    st.info("**To report biological variation data appropriately, please refer:** *Bartlett WA, Sandberg S, Carobene A, Fernandez-Calle P, Diaz-Garzon J, Coskun A, Jonker N, Galior K, Gonzales-Lao E, Moreno-Parro I, Sufrate-Vergara B, Webster C, Itkonen O, Marques-García F, Aarsand AK. A standard to report biological variation data studies - based on an expert opinion. Clin Chem Lab Med. 2024 Jul 8;63(1):52-59. doi: 10.1515/cclm-2024-0489. PMID: 38965828.*")


    st.divider()
    with st.expander("Troubleshooting", expanded=False):
        st.markdown(
            """
- **“Missing required columns”** → check mapping or rename columns then remap.
- **“'Sample/Replicate' must be numeric”** → fix non-numeric entries before mapping.
- **“All data removed during QC”** → relax switches, or untick **Enforce balanced crossed design**.
            """
        )
    st.caption("Tip: You can also download the template from the sidebar.")










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
    ready = manual_df.dropna().shape[0] >= 4
    if ready:
        user_df = _strip_ghost_index(manual_df.copy())
    # do not touch user_df otherwise; just persist edits:
    st.session_state.manual_df = manual_df


# ---- Dataset identity & remap guard ----------------------------------------
def _dataset_origin(up_file, df):
    if up_file is not None:
        # include filename to detect replacement with a different upload
        return f"upload::{up_file.name}"
    if df is not None:
        return "manual"
    return None

curr_origin = _dataset_origin(up_file, user_df)
curr_cols   = tuple(map(str, user_df.columns)) if user_df is not None else None
curr_rows   = int(len(user_df)) if user_df is not None else None

prev_origin = st.session_state.get("data_origin")
prev_cols   = st.session_state.get("data_cols")
prev_rows   = st.session_state.get("data_rows")

# Invalidate saved mapping if:
#  • source changed (upload ↔ manual, or different file name)
#  • columns changed (headers remapped)
#  • dataset cleared (user removed file → origin becomes None)
#  • row count changed (likely a different upload with same headers)
if (curr_origin != prev_origin) or (curr_cols != prev_cols) or (curr_rows != prev_rows):
    st.session_state["mapping_ok"] = False
    st.session_state["mapped_df"] = None
    st.session_state["mapped_source_cols"] = list(user_df.columns) if user_df is not None else None
    st.session_state["result_unit"] = ""   # NEW: clear unit when data fingerprint changes

# Persist current fingerprint for next run
st.session_state["data_origin"] = curr_origin
st.session_state["data_cols"]   = curr_cols
st.session_state["data_rows"]   = curr_rows
# ---------------------------------------------------------------------------

# NEW: persistent mapping state gate
if "mapping_ok" not in st.session_state:
    st.session_state["mapping_ok"] = False
if "mapped_df" not in st.session_state:
    st.session_state["mapped_df"] = None

# ——————————————————————————————————————————————————————————————
# 7. Preview & calculate button
# ——————————————————————————————————————————————————————————————
if user_df is not None:
    st.subheader("Data preview")
    st.dataframe(user_df, use_container_width=True, hide_index=True)

    if st.session_state.get("mapping_ok", False):
        st.success("A valid mapping is saved for this dataset.")
    else:
        st.info("Please map the four required columns below.")

    # — Column mapping (collapsed until user opens it) —
    # — Column mapping (open until valid) —
    with st.expander("⇢  Map / confirm columns", expanded=not st.session_state.get("mapping_ok", False)):
        st.markdown("""
    **Why mapping?** The calculator needs to know which column is the **subject**, which is the **time-point**, the **replicate index**, and the numeric **result**.  
    This ensures the ANOVA reads the hierarchy correctly (Subject → Sample/visit → Replicate).

    **How to choose**
    - **Subject identifier:** a stable participant ID (e.g., `P01`, `MRN1234`).  
    - **Sample / time-point:** the order of collections per subject (e.g., `1, 2, 3…`).  
    - **Replicate:** within-sample repeat index (`1, 2, 3…`) from the same run/tube.  
    - **Result / value:** the analyte measurement as a **single numeric** value in one unit.
    """)

        cols = list(user_df.columns)
        with st.form("map_form"):
            c1, c2 = st.columns(2)
            with c1:
                subj_sel = st.selectbox("Subject identifier",  cols,
                                        index=_default_idx("Subject", cols),
                                        help="Unique person/participant code.")
                samp_sel = st.selectbox("Sample / time-point", cols,
                                        index=_default_idx("Sample", cols),
                                        help="Visit/order per subject (1, 2, 3…).")
                result_unit = st.text_input(
                                            "Unit for Result values (e.g., mg/dL, mmol/L)",
                                            value=st.session_state.get("result_unit", ""),
                                            placeholder="e.g., mg/dL",
                                            help="This unit will be shown next to the mean concentration."
                                        )
            with c2:
                repl_sel = st.selectbox("Replicate",          cols,
                                        index=_default_idx("Replicate", cols),
                                        help="Within-sample repeat index (1, 2…).")
                res_sel  = st.selectbox("Result / value",     cols,
                                        index=_default_idx("Result", cols),
                                        help="Numeric measurement; one unit.")
                measurand_name = st.text_input(
                                            "Measurand (e.g. Glucose)",
                                            value=st.session_state.get("measurand_name", ""),
                                            placeholder="e.g. Glucose",
                                            help="Enter the name of the measurand for plot labels."
                                        )
            # Single save button (no reset button)
            confirmed = st.form_submit_button("Save mapping", use_container_width=True)
            st.caption("To change mapping later, adjust the selections above and click **Save mapping** again.")


        if confirmed:
            mapped_df, errs, warns, dropped_indices = _validate_and_build_mapping(
                user_df, subj_sel, samp_sel, repl_sel, res_sel
            )
            if errs:
                for e in errs:
                    st.error(e)
                st.session_state["mapping_ok"] = False
                st.session_state["mapped_df"] = None
                st.session_state["mapping_warnings"] = []
                st.session_state["dropped_row_indices"] = []
            else:
                if warns:
                    for w in warns:
                        st.warning(w)
                st.session_state["mapping_ok"] = True
                st.session_state["mapped_df"] = mapped_df
                st.session_state["mapping_warnings"] = warns  # Store for preprocessing log
                st.session_state["dropped_row_indices"] = dropped_indices  # Store dropped indices for gender split
                st.session_state["result_unit"] = (result_unit or "").strip()   # NEW
                st.session_state["measurand_name"] = (measurand_name or "").strip()   # NEW
                st.session_state["mapped_source_cols"] = list(user_df.columns)  # NEW
                st.success("Mapping saved ✓ — columns validated and ready to calculate.")
                st.caption(f"Preview of mapped columns (first 8 rows):")
                st.dataframe(mapped_df.head(8), use_container_width=True, hide_index=True)



    # --- BIVAC QI 1–6 manual entry (appears before Calculate) ---
    with st.expander("⇢ BIVAC – Manual entry for QI 1–6 (study design & methods)", expanded=False):
        if "qi_manual" not in st.session_state:
            st.session_state.qi_manual = {}

        for label, cfg in QI15_CHOICES.items():
            col1, col2 = st.columns([0.55, 0.45])
            with col1:
                st.markdown(f"**{label}**")
                st.caption(cfg["help"])
            with col2:
                choice = st.selectbox(
                    "Select grade",
                    options=list(cfg["options"].keys()),
                    index=0,
                    format_func=lambda k, cfg=cfg: f"{k} – {cfg['options'][k]}",
                    key=f"qi16_{label}",
                )
            st.session_state.qi_manual[label] = choice

        st.info("These selections will be included in the checklist and overall BIVAC grade.")
        st.info("**When using BIVAC, please cite the following reference:** *Aarsand AK, Røraas T, Fernandez-Calle P, Ricos C, Díaz-Garzón J, Jonker N, Perich C, González-Lao E, Carobene A, Minchinela J, Coşkun A, Simón M, Álvarez V, Bartlett WA, Fernández-Fernández P, Boned B, Braga F, Corte Z, Aslan B, Sandberg S; European Federation of Clinical Chemistry and Laboratory Medicine Working Group on Biological Variation and Task and Finish Group for the Biological Variation Database. The Biological Variation Data Critical Appraisal Checklist: A Standard for Evaluating Studies on Biological Variation. Clin Chem. 2018 Mar;64(3):501-514. doi: 10.1373/clinchem.2017.281808. Epub 2017 Dec 8. PMID: 29222339.*")


    # NEW – user option: CV-ANOVA
    estimate_cv_anova = st.checkbox("Estimate CVI with CV-ANOVA")
    st.session_state["use_cv_anova"] = estimate_cv_anova

    # ─────────────────────────────────────────────────────────────
    # Replicate-level outlier detection and user selection UI
    # ─────────────────────────────────────────────────────────────
    # Define mapped_df at this scope so ALL outlier UI sections can access it
    mapped_df = st.session_state.get("mapped_df")
    
    if st.session_state.get("mapping_ok", False) and st.session_state.get("preproc_flags", {}).get("rep_cochran", True) and mapped_df is not None:
        if mapped_df is not None and not mapped_df.empty:
            # Detect replicate-level outliers
            detected_outliers = _detect_replicate_outliers(mapped_df.copy())
            
            # Initialize selection storage if not present
            if "rep_outlier_selections" not in st.session_state:
                st.session_state["rep_outlier_selections"] = {}
            
            if detected_outliers:
                # In AUTO mode, set all selections to "remove_all" without showing UI
                if st.session_state.get("outlier_mode", "auto") == "auto":
                    for outlier in detected_outliers:
                        key = (outlier["Subject"], outlier["Sample"])
                        st.session_state["rep_outlier_selections"][key] = "remove_all"
                else:
                    # MANUAL mode: show UI for user selection
                    st.markdown("---")
                    st.markdown("#### 🔍 Replicate-level Outliers (Cochran QI 8a)")
                    st.caption(f"Found **{len(detected_outliers)}** sample(s) with unusually high replicate variance. Choose an action for each:")
                
                    for i, outlier in enumerate(detected_outliers):
                        subj = outlier["Subject"]
                        samp = outlier["Sample"]
                        variance = outlier["variance"]
                        reps = outlier["replicates"]
                        key = (subj, samp)
                        
                        # Compact card-style container
                        with st.container():
                            col1, col2 = st.columns([1, 2])
                            
                            with col1:
                                st.markdown(f"**Subject {subj}, Sample {samp}**")
                                st.caption(f"s² = {variance:.4f} | G = {outlier['G']:.3f}")
                                # Compact replicate display
                                rep_str = " · ".join([f"R{r['Replicate']}={r['Result']:.2f}" for r in reps])
                                st.code(rep_str, language=None)
                            
                            with col2:
                                # Build compact options
                                options = ["🗑️ Remove sample"]
                                for r in reps:
                                    options.append(f"✅ Keep R{r['Replicate']}")
                                options.append("⏩ Ignore")
                                
                                # Get stored selection or default
                                stored_idx = 0
                                if key in st.session_state.get("rep_outlier_selections", {}):
                                    stored_action = st.session_state["rep_outlier_selections"][key]
                                    if stored_action == "remove_all":
                                        stored_idx = 0
                                    elif stored_action == "ignore":
                                        stored_idx = len(options) - 1
                                    elif stored_action.startswith("keep_rep_"):
                                        rep_id = stored_action.replace("keep_rep_", "")
                                        for j, r in enumerate(reps):
                                            if str(r["Replicate"]) == rep_id:
                                                stored_idx = j + 1
                                                break
                                
                                selected = st.radio(
                                    "Action:",
                                    options,
                                    index=stored_idx,
                                    key=f"rep_outlier_action_{i}_{subj}_{samp}",
                                    horizontal=True,
                                    label_visibility="collapsed"
                                )
                                
                                # Map selection to action string
                                if selected.startswith("🗑️"):
                                    action = "remove_all"
                                elif selected.startswith("⏩"):
                                    action = "ignore"
                                else:
                                    # Extract replicate number from selection
                                    for r in reps:
                                        if f"R{r['Replicate']}" in selected:
                                            action = f"keep_rep_{r['Replicate']}"
                                            break
                                
                                st.session_state["rep_outlier_selections"][key] = action
                        
                        if i < len(detected_outliers) - 1:
                            st.markdown("<hr style='margin: 0.5rem 0; opacity: 0.3;'>", unsafe_allow_html=True)
            else:
                # Clear any stale selections if no outliers detected
                st.session_state["rep_outlier_selections"] = {}
        else:
            st.session_state["rep_outlier_selections"] = {}
    else:
        # Clear selections if Cochran is disabled
        st.session_state["rep_outlier_selections"] = {}

    # ─────────────────────────────────────────────────────────────
    # Bartlett Outlier Detection & Selection UI (QI 8a - Bartlett)
    # ─────────────────────────────────────────────────────────────
    if preproc_flags.get("rep_bartlett", False) and st.session_state.get("mapping_ok", False) and mapped_df is not None:
        # Initialize session state for Bartlett selections
        if "rep_bartlett_selections" not in st.session_state:
            st.session_state["rep_bartlett_selections"] = {}
        
        # Apply Cochran exclusions first (Bartlett should run on Cochran-cleaned data)
        df_for_bartlett = mapped_df.copy()
        cochran_selections = st.session_state.get("rep_outlier_selections", {})
        for (subj, samp), action in cochran_selections.items():
            if action == "remove_all":
                df_for_bartlett = df_for_bartlett[~((df_for_bartlett["Subject"] == subj) & (df_for_bartlett["Sample"] == samp))]
            elif action.startswith("keep_rep_"):
                rep_to_keep = action.replace("keep_rep_", "")
                try:
                    rep_to_keep = int(float(rep_to_keep))
                except ValueError:
                    pass
                mask = (df_for_bartlett["Subject"] == subj) & (df_for_bartlett["Sample"] == samp)
                keep_mask = mask & (df_for_bartlett["Replicate"] == rep_to_keep)
                df_for_bartlett = df_for_bartlett[~mask | keep_mask]
            # "ignore" action: keep data as-is
        
        # Detect Bartlett outliers on Cochran-cleaned data
        bartlett_outliers = _detect_bartlett_outliers(df_for_bartlett)
        
        if bartlett_outliers:
            if st.session_state.get("outlier_mode", "auto") == "auto":
                for outlier in bartlett_outliers:
                    key = (outlier["Subject"], outlier["Sample"])
                    st.session_state["rep_bartlett_selections"][key] = "remove_all"
            else:
                st.markdown("---")
                st.markdown("#### 🔬 Bartlett Iterative Variance Outlier Exclusion (QI 8a)")
                st.info(
                    f"**Iterative Process:** Found **{len(bartlett_outliers)}** sample(s) to review. "
                    "Bartlett's test removes the **highest-variance sample**, re-tests, and repeats until homogeneity is achieved. "
                    "The order below reflects the removal sequence."
                )
            
                for i, outlier in enumerate(bartlett_outliers):
                    subj = outlier["Subject"]
                    samp = outlier["Sample"]
                    variance = outlier["variance"]
                    reps = outlier["replicates"]
                    key = (subj, samp)
                    
                    # Step indicator with visual hierarchy
                    step_num = i + 1
                    with st.container():
                    # Step header
                        st.markdown(
                            f"<div style='background: linear-gradient(90deg, #667eea22 0%, transparent 100%); "
                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #667eea; margin-bottom: 0.5rem;'>"
                            f"<b>Step {step_num}</b> — Subject {subj}, Sample {samp}</div>",
                            unsafe_allow_html=True
                        )
                        
                        col1, col2 = st.columns([1, 2])
                        
                        with col1:
                            st.caption(f"**Variance:** s² = {variance:.4f}")
                            st.caption(f"**Bartlett p:** {outlier['p']:.4f}")
                            # Compact replicate display
                            rep_str = " · ".join([f"R{r['Replicate']}={r['Result']:.2f}" for r in reps])
                            st.code(rep_str, language=None)
                        
                        with col2:
                            # Build compact options with step-aware labels
                            options = [f"🗑️ Remove (proceed to Step {step_num + 1})"]
                            for r in reps:
                                options.append(f"✅ Keep only R{r['Replicate']}")
                            options.append("⏩ Ignore (accept heterogeneity)")
                            
                            # Get stored selection or default
                            stored_idx = 0
                            if key in st.session_state.get("rep_bartlett_selections", {}):
                                stored_action = st.session_state["rep_bartlett_selections"][key]
                                if stored_action == "remove_all":
                                    stored_idx = 0
                                elif stored_action == "ignore":
                                    stored_idx = len(options) - 1
                                elif stored_action.startswith("keep_rep_"):
                                    rep_id = stored_action.replace("keep_rep_", "")
                                    for j, r in enumerate(reps):
                                        if str(r["Replicate"]) == rep_id:
                                            stored_idx = j + 1
                                            break
                            
                            selected = st.radio(
                                "Action:",
                                options,
                                index=stored_idx,
                                key=f"bartlett_outlier_action_{i}_{subj}_{samp}",
                                horizontal=True,
                                label_visibility="collapsed"
                            )
                            
                            # Map selection to action string
                            if selected.startswith("🗑️"):
                                action = "remove_all"
                            elif selected.startswith("⏩"):
                                action = "ignore"
                            else:
                                # Extract replicate number from selection
                                for r in reps:
                                    if f"R{r['Replicate']}" in selected:
                                        action = f"keep_rep_{r['Replicate']}"
                                        break
                            
                            st.session_state["rep_bartlett_selections"][key] = action
                
                # Arrow connector between steps
                if i < len(bartlett_outliers) - 1:
                    st.markdown(
                        "<div style='text-align: center; color: #667eea; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>",
                        unsafe_allow_html=True
                    )
            
                # Final state indicator (outside loop)
                st.markdown(
                    "<div style='background: #22c55e22; padding: 0.5rem; border-left: 3px solid #22c55e; margin-top: 0.5rem;'>"
                    "✅ <b>Homogeneity Target</b> — After processing all steps above, variances should be homogeneous.</div>",
                    unsafe_allow_html=True
                )
        else:
                    # Clear any stale selections if no outliers detected
            st.session_state["rep_bartlett_selections"] = {}
    else:
        # Clear selections if Bartlett exclusion is disabled
        st.session_state["rep_bartlett_selections"] = {}


    # ─────────────────────────────────────────────────────────────
    # QI 8b – Sample‑level Cochran Outlier Detection & Selection UI
    # ─────────────────────────────────────────────────────────────
    if preproc_flags.get("samp_cochran", True) and st.session_state.get("mapping_ok", False) and user_df is not None:
        if "samp_cochran_selections" not in st.session_state:
            st.session_state["samp_cochran_selections"] = {}
            
        # 1. Clean data based on previous steps (Replicate Cochran & Bartlett)
        df_clean_sc = mapped_df.copy()
        
        # Apply Rep Cochran
        rc_sels = st.session_state.get("rep_outlier_selections", {})
        for (subj, samp), action in rc_sels.items():
            if action == "remove_all":
                df_clean_sc = df_clean_sc[~((df_clean_sc["Subject"] == subj) & (df_clean_sc["Sample"] == samp))]
            elif action.startswith("keep_rep_"):
                try:
                    rep_to_keep = int(float(action.replace("keep_rep_", "")))
                    mask = (df_clean_sc["Subject"] == subj) & (df_clean_sc["Sample"] == samp)
                    keep_mask = mask & (df_clean_sc["Replicate"] == rep_to_keep)
                    df_clean_sc = df_clean_sc[~mask | keep_mask]
                except: pass

        # Apply Rep Bartlett
        rb_sels = st.session_state.get("rep_bartlett_selections", {})
        for (subj, samp), action in rb_sels.items():
            if action == "remove_all":
                df_clean_sc = df_clean_sc[~((df_clean_sc["Subject"] == subj) & (df_clean_sc["Sample"] == samp))]
            elif action.startswith("keep_rep_"):
                try:
                    rep_to_keep = int(float(action.replace("keep_rep_", "")))
                    mask = (df_clean_sc["Subject"] == subj) & (df_clean_sc["Sample"] == samp)
                    keep_mask = mask & (df_clean_sc["Replicate"] == rep_to_keep)
                    df_clean_sc = df_clean_sc[~mask | keep_mask]
                except: pass
        
        # 2. Detect Sample Outliers
        samp_outliers = _detect_sample_cochran_outliers(df_clean_sc)
        
        if samp_outliers:
            if st.session_state.get("outlier_mode", "auto") == "auto":
                for out in samp_outliers:
                    st.session_state["samp_cochran_selections"][(out["Subject"], out["Sample"])] = "remove"
            else:
                st.markdown("---")
                st.markdown("#### 🔬 Sample-level Outliers (Cochran QI 8b)")
                st.info(f"**Iterative Process:** Found **{len(samp_outliers)}** outlier sample(s).")
                
                for i, outlier in enumerate(samp_outliers):
                    subj = outlier["Subject"]
                    samp = outlier["Sample"]
                    variance = outlier["variance"]
                    key = (subj, samp)
                    step_num = i + 1
                    
                    with st.container():
                        st.markdown(
                            f"<div style='background: linear-gradient(90deg, #8b5cf622 0%, transparent 100%); "
                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #8b5cf6; margin-bottom: 0.5rem;'>"
                            f"<b>Step {step_num}</b> — Subject {subj}, Sample {samp}</div>",
                            unsafe_allow_html=True
                        )
                        col1, col2 = st.columns([1, 2])
                        with col1:
                            st.caption(f"**Variance:** s² = {variance:.4f}")
                            st.caption(f"**G:** {outlier['G']:.3f} (Crit: {outlier['G_crit']:.3f})")
                        with col2:
                            options = ["🗑️ Remove Sample", "⏩ Ignore (Stop)"]
                            stored_idx = 0
                            if key in st.session_state.get("samp_cochran_selections", {}):
                                if st.session_state["samp_cochran_selections"][key] == "ignore":
                                    stored_idx = 1
                            
                            selected = st.radio("Action:", options, index=stored_idx, 
                                              key=f"samp_coch_{i}_{subj}_{samp}", horizontal=True, label_visibility="collapsed")
                            
                            action = "remove" if selected.startswith("🗑️") else "ignore"
                            st.session_state["samp_cochran_selections"][key] = action
                    
                    if i < len(samp_outliers) - 1:
                        st.markdown("<div style='text-align: center; color: #8b5cf6; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)
                st.markdown("---")
        else:
             st.session_state["samp_cochran_selections"] = {}
    else:
        st.session_state["samp_cochran_selections"] = {}



    # ─────────────────────────────────────────────────────────────
    # QI 10 Within-subject Bartlett Outlier Detection & Selection UI
    # ─────────────────────────────────────────────────────────────
    if preproc_flags.get("wp_bartlett_exclusion", False) and st.session_state.get("mapping_ok", False) and user_df is not None:
        # Initialize session state for QI 10 selections
        if "wp_bartlett_selections" not in st.session_state:
            st.session_state["wp_bartlett_selections"] = {}
        
        # Apply previous exclusions before detecting QI 10 outliers
        # QI 10 should run on data after Cochran and QI 8a Bartlett are applied
        df_for_wp_bartlett = mapped_df.copy()
        
        # Apply Cochran selections
        cochran_sels = st.session_state.get("rep_outlier_selections", {})
        for (subj, samp), action in cochran_sels.items():
            if action == "remove_all":
                df_for_wp_bartlett = df_for_wp_bartlett[~((df_for_wp_bartlett["Subject"] == subj) & (df_for_wp_bartlett["Sample"] == samp))]
            elif action.startswith("keep_rep_"):
                rep_to_keep = action.replace("keep_rep_", "")
                try:
                    rep_to_keep = int(float(rep_to_keep))
                except ValueError:
                    pass
                mask = (df_for_wp_bartlett["Subject"] == subj) & (df_for_wp_bartlett["Sample"] == samp)
                keep_mask = mask & (df_for_wp_bartlett["Replicate"] == rep_to_keep)
                df_for_wp_bartlett = df_for_wp_bartlett[~mask | keep_mask]
        
        # Apply QI 8a Bartlett selections
        bartlett_sels = st.session_state.get("rep_bartlett_selections", {})
        for (subj, samp), action in bartlett_sels.items():
            if action == "remove_all":
                df_for_wp_bartlett = df_for_wp_bartlett[~((df_for_wp_bartlett["Subject"] == subj) & (df_for_wp_bartlett["Sample"] == samp))]
            elif action.startswith("keep_rep_"):
                rep_to_keep = action.replace("keep_rep_", "")
                try:
                    rep_to_keep = int(float(rep_to_keep))
                except ValueError:
                    pass
                mask = (df_for_wp_bartlett["Subject"] == subj) & (df_for_wp_bartlett["Sample"] == samp)
                keep_mask = mask & (df_for_wp_bartlett["Replicate"] == rep_to_keep)
                df_for_wp_bartlett = df_for_wp_bartlett[~mask | keep_mask]
                
        # Apply Sample Cochran selections
        sc_sels = st.session_state.get("samp_cochran_selections", {})
        for (subj, samp), action in sc_sels.items():
             if action == "remove":
                df_for_wp_bartlett = df_for_wp_bartlett[~((df_for_wp_bartlett["Subject"] == subj) & (df_for_wp_bartlett["Sample"] == samp))]
        
        # Detect QI 10 outliers
        wp_bartlett_outliers = _detect_wp_bartlett_outliers(df_for_wp_bartlett)
        
        if wp_bartlett_outliers:
            # In AUTO mode, set all selections to "remove" without showing UI
            if st.session_state.get("outlier_mode", "auto") == "auto":
                for outlier in wp_bartlett_outliers:
                    st.session_state["wp_bartlett_selections"][outlier["Subject"]] = "remove"
            else:
                # MANUAL mode: show UI for user selection
                st.markdown("---")
                st.markdown("#### 🧬 QI 10: Within-Subject Variance Outliers (Bartlett)")
                st.info(
                    f"**Iterative Process:** Found **{len(wp_bartlett_outliers)}** subject(s) causing within-subject variance heterogeneity. "
                    "The order below reflects the removal sequence."
                )
            
                for i, outlier in enumerate(wp_bartlett_outliers):
                    subj = outlier["Subject"]
                    variance = outlier["within_subject_variance"]
                    sample_means = outlier["sample_means"]
                    step_num = i + 1
                    
                    with st.container():
                                st.markdown(
                                    f"<div style='background: linear-gradient(90deg, #f59e0b22 0%, transparent 100%); "
                                    f"padding: 0.3rem 0.5rem; border-left: 3px solid #f59e0b; margin-bottom: 0.5rem;'>"
                                    f"<b>Step {step_num}</b> — Subject {subj}</div>",
                                    unsafe_allow_html=True
                                )
                                
                                col1, col2 = st.columns([1, 2])
                                
                                with col1:
                                    st.caption(f"**Within-subject s²:** {variance:.4f}")
                                    st.caption(f"**Bartlett p:** {outlier['p']:.4f}")
                                    # Show sample means
                                    means_str = " · ".join([f"S{sm['Sample']}={sm['Result']:.2f}" for sm in sample_means])
                                    st.code(means_str, language=None)
                                
                                with col2:
                                    options = [f"🗑️ Remove subject (proceed to Step {step_num + 1})", "⏩ Ignore (accept heterogeneity)"]
                                    
                                    stored_idx = 0
                                    if subj in st.session_state.get("wp_bartlett_selections", {}):
                                        stored_action = st.session_state["wp_bartlett_selections"][subj]
                                        if stored_action == "remove":
                                            stored_idx = 0
                                        elif stored_action == "ignore":
                                            stored_idx = 1
                                    
                                    selected = st.radio(
                                        "Action:",
                                        options,
                                        index=stored_idx,
                                        key=f"wp_bartlett_action_{i}_{subj}",
                                        horizontal=True,
                                        label_visibility="collapsed"
                                    )
                                    
                                    if selected.startswith("🗑️"):
                                        action = "remove"
                                    else:
                                        action = "ignore"
                                    
                                    st.session_state["wp_bartlett_selections"][subj] = action
                            
                    if i < len(wp_bartlett_outliers) - 1:
                        st.markdown(
                            "<div style='text-align: center; color: #f59e0b; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>",
                            unsafe_allow_html=True
                        )
                
                st.markdown(
                    "<div style='background: #22c55e22; padding: 0.5rem; border-left: 3px solid #22c55e; margin-top: 0.5rem;'>"
                    "✅ <b>Homogeneity Target</b> — After processing all steps, within-subject variances should be homogeneous.</div>",
                    unsafe_allow_html=True
                )
        else:
                    st.session_state["wp_bartlett_selections"] = {}
 
    # ─────────────────────────────────────────────────────────────
    # QI 8c – Reed Between-Subject Outlier Detection & Selection UI
    # ─────────────────────────────────────────────────────────────
    if preproc_flags.get("reed", True) and st.session_state.get("mapping_ok", False) and user_df is not None:
        if "reed_selections" not in st.session_state:
            st.session_state["reed_selections"] = {}
            
        # 1. Clean data based on ALL previous steps (8a -> 8b -> 10)
        df_clean_reed = mapped_df.copy()
        
        # Apply Rep Cochran
        rc_sels = st.session_state.get("rep_outlier_selections", {})
        for (subj, samp), action in rc_sels.items():
            if action == "remove_all":
                df_clean_reed = df_clean_reed[~((df_clean_reed["Subject"] == subj) & (df_clean_reed["Sample"] == samp))]
            elif action.startswith("keep_rep_"):
                try:
                    rep_to_keep = int(float(action.replace("keep_rep_", "")))
                    mask = (df_clean_reed["Subject"] == subj) & (df_clean_reed["Sample"] == samp)
                    keep_mask = mask & (df_clean_reed["Replicate"] == rep_to_keep)
                    df_clean_reed = df_clean_reed[~mask | keep_mask]
                except: pass

        # Apply Rep Bartlett
        rb_sels = st.session_state.get("rep_bartlett_selections", {})
        for (subj, samp), action in rb_sels.items():
            if action == "remove_all":
                df_clean_reed = df_clean_reed[~((df_clean_reed["Subject"] == subj) & (df_clean_reed["Sample"] == samp))]
            elif action.startswith("keep_rep_"):
                try:
                    rep_to_keep = int(float(action.replace("keep_rep_", "")))
                    mask = (df_clean_reed["Subject"] == subj) & (df_clean_reed["Sample"] == samp)
                    keep_mask = mask & (df_clean_reed["Replicate"] == rep_to_keep)
                    df_clean_reed = df_clean_reed[~mask | keep_mask]
                except: pass
                
        # Apply Sample Cochran
        sc_sels = st.session_state.get("samp_cochran_selections", {})
        for (subj, samp), action in sc_sels.items():
             if action == "remove":
                df_clean_reed = df_clean_reed[~((df_clean_reed["Subject"] == subj) & (df_clean_reed["Sample"] == samp))]

        # Apply QI 10 Bartlett (New!)
        wb_sels = st.session_state.get("wp_bartlett_selections", {})
        for subj, action in wb_sels.items():
            if action == "remove":
                df_clean_reed = df_clean_reed[df_clean_reed["Subject"] != subj]

        # 2. Detect Reed Outliers
        reed_outliers = _detect_reed_outliers(df_clean_reed)
        
        if reed_outliers:
            if st.session_state.get("outlier_mode", "auto") == "auto":
                for out in reed_outliers:
                    st.session_state["reed_selections"][out["Subject"]] = "remove"
            else:
                st.markdown("---")
                st.markdown("#### 🔬 Reed Outliers (QI 8c)")
                st.info(f"**Between-subject check:** Found **{len(reed_outliers)}** outlier subject(s) with extreme means.")
                
                for i, outlier in enumerate(reed_outliers):
                    subj = outlier["Subject"]
                    mean_val = outlier["mean"]
                    reason = outlier["reason"]
                    step_num = i + 1
                    
                    with st.container():
                        st.markdown(
                            f"<div style='background: linear-gradient(90deg, #ef444422 0%, transparent 100%); "
                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #ef4444; margin-bottom: 0.5rem;'>"
                            f"<b>Step {step_num}</b> — Subject {subj}</div>",
                            unsafe_allow_html=True
                        )
                        col1, col2 = st.columns([1, 2])
                        with col1:
                            st.caption(f"**Mean:** {mean_val:.4f}")
                            st.caption(f"**Reason:** {reason}")
                        with col2:
                            options = ["🗑️ Remove Subject", "⏩ Ignore (Stop)"]
                            stored_idx = 0
                            if subj in st.session_state.get("reed_selections", {}):
                                if st.session_state["reed_selections"][subj] == "ignore":
                                    stored_idx = 1
                            
                            selected = st.radio("Action:", options, index=stored_idx, 
                                              key=f"reed_{i}_{subj}", horizontal=True, label_visibility="collapsed")
                            
                            action = "remove" if selected.startswith("🗑️") else "ignore"
                            st.session_state["reed_selections"][subj] = action
                    
                    if i < len(reed_outliers) - 1:
                        st.markdown("<div style='text-align: center; color: #ef4444; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)
                st.markdown("---")
        else:
            if st.session_state.get("outlier_mode", "auto") != "auto":
                st.markdown("#### 🔬 Reed Outliers (QI 8c)")
                st.success("✅ No outliers detected (or masked by prior exclusion steps).")
                st.markdown("---")
            st.session_state["reed_selections"] = {}

    # ─────────────────────────────────────────────────────────────
    # QI 7 – Steady-State Drift Detection & Selection UI
    # ─────────────────────────────────────────────────────────────
    if preproc_flags.get("drift", True) and st.session_state.get("mapping_ok", False) and user_df is not None:
        if "drift_selections" not in st.session_state:
            st.session_state["drift_selections"] = {}
            
        # 1. Clean data based on ALL previous steps
        # Cascade: Rep Cochran -> Rep Bartlett -> Samp Cochran -> WP Bartlett -> Reed
        df_clean_drift = mapped_df.copy()
        
        # Apply Rep Cochran
        rc_sels = st.session_state.get("rep_outlier_selections", {})
        for (subj, samp), action in rc_sels.items():
            if action == "remove_all":
                df_clean_drift = df_clean_drift[~((df_clean_drift["Subject"] == subj) & (df_clean_drift["Sample"] == samp))]
            elif action.startswith("keep_rep_"):
                try:
                    rep_to_keep = int(float(action.replace("keep_rep_", "")))
                    mask = (df_clean_drift["Subject"] == subj) & (df_clean_drift["Sample"] == samp)
                    keep_mask = mask & (df_clean_drift["Replicate"] == rep_to_keep)
                    df_clean_drift = df_clean_drift[~mask | keep_mask]
                except: pass

        # Apply Rep Bartlett
        rb_sels = st.session_state.get("rep_bartlett_selections", {})
        for (subj, samp), action in rb_sels.items():
            if action == "remove_all":
                df_clean_drift = df_clean_drift[~((df_clean_drift["Subject"] == subj) & (df_clean_drift["Sample"] == samp))]
            elif action.startswith("keep_rep_"):
                try:
                    rep_to_keep = int(float(action.replace("keep_rep_", "")))
                    mask = (df_clean_drift["Subject"] == subj) & (df_clean_drift["Sample"] == samp)
                    keep_mask = mask & (df_clean_drift["Replicate"] == rep_to_keep)
                    df_clean_drift = df_clean_drift[~mask | keep_mask]
                except: pass
                
        # Apply Sample Cochran
        sc_sels = st.session_state.get("samp_cochran_selections", {})
        for (subj, samp), action in sc_sels.items():
             if action == "remove":
                df_clean_drift = df_clean_drift[~((df_clean_drift["Subject"] == subj) & (df_clean_drift["Sample"] == samp))]

        # Apply QI 10 Bartlett
        wb_sels = st.session_state.get("wp_bartlett_selections", {})
        for subj, action in wb_sels.items():
            if action == "remove":
                df_clean_drift = df_clean_drift[df_clean_drift["Subject"] != subj]

        # Apply Reed
        reed_sels = st.session_state.get("reed_selections", {})
        for subj, action in reed_sels.items():
            if action == "remove":
                df_clean_drift = df_clean_drift[df_clean_drift["Subject"] != subj]

        # 2. Detect Drift Outliers
        drift_outliers = _detect_drift_outliers(df_clean_drift)
        
        if drift_outliers:
            if st.session_state.get("outlier_mode", "auto") == "auto":
                for out in drift_outliers:
                    st.session_state["drift_selections"][out["Subject"]] = "remove"
            else:
                st.markdown("---")
                st.markdown("#### 📉 Steady-State Drift (QI 7)")
                st.info(f"**Trend check:** Found **{len(drift_outliers)}** subject(s) with significant drift over time.")
                
                for i, outlier in enumerate(drift_outliers):
                    subj = outlier["Subject"]
                    slope = outlier["slope"]
                    p_val = outlier["p"]
                    step_num = i + 1
                    
                    with st.container():
                        st.markdown(
                            f"<div style='background: linear-gradient(90deg, #8b5cf622 0%, transparent 100%); "
                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #8b5cf6; margin-bottom: 0.5rem;'>"
                            f"<b>Step {step_num}</b> — Subject {subj}</div>",
                            unsafe_allow_html=True
                        )
                        col1, col2 = st.columns([1, 2])
                        with col1:
                            st.caption(f"**Slope:** {slope:.4f}")
                            st.caption(f"**p-value:** {p_val:.4f}")
                        with col2:
                            options = ["🗑️ Remove Subject", "⏩ Ignore (Accept Trend)"]
                            stored_idx = 0
                            if subj in st.session_state.get("drift_selections", {}):
                                if st.session_state["drift_selections"][subj] == "ignore":
                                    stored_idx = 1
                            
                            selected = st.radio("Action:", options, index=stored_idx, 
                                              key=f"drift_{i}_{subj}", horizontal=True, label_visibility="collapsed")
                            
                            action = "remove" if selected.startswith("🗑️") else "ignore"
                            st.session_state["drift_selections"][subj] = action
                    
                    if i < len(drift_outliers) - 1:
                        st.markdown("<div style='text-align: center; color: #8b5cf6; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)
                st.markdown("---")
        else:
             st.session_state["drift_selections"] = {}
    else:
        st.session_state["drift_selections"] = {}

    # ─────────────────────────────────────────────────────────────
    # Gender-based calculations UI
    # ─────────────────────────────────────────────────────────────
    gender_based = st.checkbox("Gender based calculations", value=False)
    st.session_state["gender_based"] = gender_based

    gender_settings_ok = True
    if gender_based:
        st.markdown("### Gender mapping")

        cols_all = list(user_df.columns)  # raw dataset columns

        gender_mode = st.radio(
            "How is gender encoded in your dataset?",
            ["Single gender column", "Two indicator columns (Male/Female)"],
            horizontal=True
        )
        st.session_state["gender_mode"] = gender_mode

        if gender_mode == "Single gender column":
            gender_col = st.selectbox("Gender column", cols_all, index=0)
            # Pull unique values for easy mapping
            uniq = pd.Series(user_df[gender_col]).dropna().unique().tolist()
            if len(uniq) == 0:
                st.error("Selected gender column has no values.")
                gender_settings_ok = False
            else:
                male_value = st.selectbox("Which variable means **Male**?", uniq, index=0)
                # Try to pick a different default for female if possible
                f_index = 1 if len(uniq) > 1 else 0
                female_value = st.selectbox("Which variable means **Female**?", uniq, index=f_index)

                if male_value == female_value:
                    st.error("Male and Female values must be different.")
                    gender_settings_ok = False

                st.session_state["gender_col"] = gender_col
                st.session_state["male_value"] = male_value
                st.session_state["female_value"] = female_value

        else:
            male_col = st.selectbox("Male indicator column", cols_all, index=0)
            female_col = st.selectbox("Female indicator column", cols_all, index=min(1, len(cols_all)-1))

            if male_col == female_col:
                st.error("Male and Female indicator columns must be different.")
                gender_settings_ok = False

            st.session_state["male_col"] = male_col
            st.session_state["female_col"] = female_col

        st.caption(
            "Rows that cannot be mapped unambiguously to Male/Female (e.g., missing/unknown/both indicators) "
            "will be excluded **only for gender-split analysis**."
        )

    st.session_state["gender_settings_ok"] = gender_settings_ok

    # ─────────────────────────────────────────────────────────────
    # Gender-specific replicate outlier detection UI
    # ─────────────────────────────────────────────────────────────
    
    # Helper function to apply Cochran/Bartlett selections to a dataframe
    # Defined at this scope so all gender-based outlier sections can use it
    def _apply_cochran_selections(df, selections):
        for (subj, samp), action in selections.items():
            if action == "remove_all" or action == "remove":
                df = df[~((df["Subject"] == subj) & (df["Sample"] == samp))]
            elif action.startswith("keep_rep_"):
                rep_to_keep = action.replace("keep_rep_", "")
                try:
                    rep_to_keep = int(float(rep_to_keep))
                except ValueError:
                    pass
                mask = (df["Subject"] == subj) & (df["Sample"] == samp)
                keep_mask = mask & (df["Replicate"] == rep_to_keep)
                df = df[~mask | keep_mask]
        return df
    
    if (gender_based and gender_settings_ok and
        st.session_state.get("mapping_ok", False) and
        st.session_state.get("preproc_flags", {}).get("rep_cochran", True)):

        mapped_df = st.session_state.get("mapped_df")
        if mapped_df is not None and not mapped_df.empty:
            # Build gender-split datasets
            df_with_gender = mapped_df.copy()
            mode = st.session_state.get("gender_mode", "Single gender column")

            if mode == "Single gender column":
                gcol = st.session_state.get("gender_col")
                mv = st.session_state.get("male_value")
                fv = st.session_state.get("female_value")
                if gcol and gcol in user_df.columns:
                    df_with_gender["Gender"] = _derive_gender_single_column(user_df, gcol, mv, fv)
            else:
                mcol = st.session_state.get("male_col")
                fcol = st.session_state.get("female_col")
                if mcol and fcol and mcol in user_df.columns and fcol in user_df.columns:
                    df_with_gender["Gender"] = _derive_gender_two_indicators(user_df, mcol, fcol)

            if "Gender" in df_with_gender.columns:
                male_df = df_with_gender[df_with_gender["Gender"] == "Male"].drop(columns=["Gender"]).copy()
                female_df = df_with_gender[df_with_gender["Gender"] == "Female"].drop(columns=["Gender"]).copy()

                # Initialize ALL gender-specific selection storage
                if "rep_outlier_selections_male" not in st.session_state:
                    st.session_state["rep_outlier_selections_male"] = {}
                if "rep_outlier_selections_female" not in st.session_state:
                    st.session_state["rep_outlier_selections_female"] = {}
                if "rep_bartlett_selections_male" not in st.session_state:
                    st.session_state["rep_bartlett_selections_male"] = {}
                if "rep_bartlett_selections_female" not in st.session_state:
                    st.session_state["rep_bartlett_selections_female"] = {}
                if "samp_cochran_selections_male" not in st.session_state:
                    st.session_state["samp_cochran_selections_male"] = {}
                if "samp_cochran_selections_female" not in st.session_state:
                    st.session_state["samp_cochran_selections_female"] = {}
                if "wp_bartlett_selections_male" not in st.session_state:
                    st.session_state["wp_bartlett_selections_male"] = {}
                if "wp_bartlett_selections_female" not in st.session_state:
                    st.session_state["wp_bartlett_selections_female"] = {}
                if "reed_selections_male" not in st.session_state:
                    st.session_state["reed_selections_male"] = {}
                if "reed_selections_female" not in st.session_state:
                    st.session_state["reed_selections_female"] = {}
                if "drift_selections_male" not in st.session_state:
                    st.session_state["drift_selections_male"] = {}
                if "drift_selections_female" not in st.session_state:
                    st.session_state["drift_selections_female"] = {}

                # ═══════════════════════════════════════════════════════════════
                # FEMALE - All Outlier Steps
                # ═══════════════════════════════════════════════════════════════
                if not female_df.empty:
                    if st.session_state.get("outlier_mode", "auto") != "auto":
                        st.markdown("---")
                        st.markdown("## 👩 Female Outlier Analysis")
                        st.markdown("---")

                    # ─────────────────────────────────────────────────────────────
                    # Female: Replicate-level Cochran (QI 8a)
                    # ─────────────────────────────────────────────────────────────
                    female_outliers = _detect_replicate_outliers(female_df)
                    if female_outliers:
                        if st.session_state.get("outlier_mode", "auto") == "auto":
                            for outlier in female_outliers:
                                key = (outlier["Subject"], outlier["Sample"])
                                st.session_state["rep_outlier_selections_female"][key] = "remove_all"
                        else:
                            st.markdown("#### 🔍 Female: Replicate-level Outliers (Cochran QI 8a)")
                            st.caption(f"Found **{len(female_outliers)}** sample(s) with high variance.")

                            for i, outlier in enumerate(female_outliers):
                                subj, samp = outlier["Subject"], outlier["Sample"]
                                variance, reps = outlier["variance"], outlier["replicates"]
                                key = (subj, samp)

                                with st.container():
                                    col1, col2 = st.columns([1, 2])

                                    with col1:
                                        st.markdown(f"**Subject {subj}, Sample {samp}**")
                                        st.caption(f"s² = {variance:.4f} | G = {outlier['G']:.3f}")
                                        rep_str = " · ".join([f"R{r['Replicate']}={r['Result']:.2f}" for r in reps])
                                        st.code(rep_str, language=None)

                                    with col2:
                                        options = ["🗑️ Remove sample"]
                                        for r in reps:
                                            options.append(f"✅ Keep R{r['Replicate']}")
                                        options.append("⏩ Ignore")

                                        stored_idx = 0
                                        if key in st.session_state.get("rep_outlier_selections_female", {}):
                                            stored_action = st.session_state["rep_outlier_selections_female"][key]
                                            if stored_action == "remove_all": stored_idx = 0
                                            elif stored_action == "ignore": stored_idx = len(options) - 1
                                            elif stored_action.startswith("keep_rep_"):
                                                for j, r in enumerate(reps):
                                                    if str(r["Replicate"]) == stored_action.replace("keep_rep_", ""):
                                                        stored_idx = j + 1; break

                                        selected = st.radio("Action:", options, index=stored_idx,
                                            key=f"female_outlier_{i}_{subj}_{samp}", horizontal=True, label_visibility="collapsed")

                                        if selected.startswith("🗑️"): action = "remove_all"
                                        elif selected.startswith("⏩"): action = "ignore"
                                        else:
                                            for r in reps:
                                                if f"R{r['Replicate']}" in selected:
                                                    action = f"keep_rep_{r['Replicate']}"; break
                                        st.session_state["rep_outlier_selections_female"][key] = action

                                if i < len(female_outliers) - 1:
                                    st.markdown("<hr style='margin: 0.5rem 0; opacity: 0.3;'>", unsafe_allow_html=True)
                            st.markdown("---")
                    else:
                        st.session_state["rep_outlier_selections_female"] = {}

                    # ─────────────────────────────────────────────────────────────
                    # Female: Replicate-level Bartlett (QI 8a)
                    # ─────────────────────────────────────────────────────────────
                    if preproc_flags.get("rep_bartlett_exclusion", False):
                        female_df_bartlett = _apply_cochran_selections(female_df.copy(), st.session_state.get("rep_outlier_selections_female", {}))
                        female_bartlett_outliers = _detect_bartlett_outliers(female_df_bartlett)

                        if female_bartlett_outliers:
                            if st.session_state.get("outlier_mode", "auto") == "auto":
                                for outlier in female_bartlett_outliers:
                                    key = (outlier["Subject"], outlier["Sample"])
                                    st.session_state["rep_bartlett_selections_female"][key] = "remove_all"
                            else:
                                st.markdown("#### 🔬 Female: Bartlett Iterative Outliers (QI 8a)")
                                st.info(f"**After Cochran cleanup:** Found **{len(female_bartlett_outliers)}** sample(s) for iterative Bartlett exclusion.")

                                for i, outlier in enumerate(female_bartlett_outliers):
                                    subj, samp = outlier["Subject"], outlier["Sample"]
                                    variance, reps = outlier["variance"], outlier["replicates"]
                                    key = (subj, samp)
                                    step_num = i + 1

                                    with st.container():
                                        st.markdown(
                                            f"<div style='background: linear-gradient(90deg, #ec489922 0%, transparent 100%); "
                                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #ec4899; margin-bottom: 0.5rem;'>"
                                            f"<b>Step {step_num}</b> — Subject {subj}, Sample {samp}</div>",
                                            unsafe_allow_html=True
                                        )
                                        col1, col2 = st.columns([1, 2])

                                        with col1:
                                            st.caption(f"**Variance:** s² = {variance:.4f}")
                                            rep_str = " · ".join([f"R{r['Replicate']}={r['Result']:.2f}" for r in reps])
                                            st.code(rep_str, language=None)

                                        with col2:
                                            options = [f"🗑️ Remove (proceed to Step {step_num + 1})"]
                                            for r in reps:
                                                options.append(f"✅ Keep only R{r['Replicate']}")
                                            options.append("⏩ Ignore (accept heterogeneity)")

                                            stored_idx = 0
                                            if key in st.session_state.get("rep_bartlett_selections_female", {}):
                                                stored_action = st.session_state["rep_bartlett_selections_female"][key]
                                                if stored_action == "remove_all": stored_idx = 0
                                                elif stored_action == "ignore": stored_idx = len(options) - 1
                                                elif stored_action.startswith("keep_rep_"):
                                                    for j, r in enumerate(reps):
                                                        if str(r["Replicate"]) == stored_action.replace("keep_rep_", ""):
                                                            stored_idx = j + 1; break

                                            selected = st.radio("Action:", options, index=stored_idx,
                                                key=f"female_bartlett_{i}_{subj}_{samp}", horizontal=True, label_visibility="collapsed")

                                            if selected.startswith("🗑️"): action = "remove_all"
                                            elif selected.startswith("⏩"): action = "ignore"
                                            else:
                                                for r in reps:
                                                    if f"R{r['Replicate']}" in selected:
                                                        action = f"keep_rep_{r['Replicate']}"; break
                                            st.session_state["rep_bartlett_selections_female"][key] = action

                                    if i < len(female_bartlett_outliers) - 1:
                                        st.markdown("<div style='text-align: center; color: #ec4899; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)
                                st.markdown(
                                    "<div style='background: #22c55e22; padding: 0.5rem; border-left: 3px solid #22c55e; margin-top: 0.5rem;'>"
                                    "✅ <b>Homogeneity Target</b> — After processing all steps above, variances should be homogeneous.</div>",
                                    unsafe_allow_html=True
                                )
                                st.markdown("---")
                        else:
                            st.session_state["rep_bartlett_selections_female"] = {}

                    # ─────────────────────────────────────────────────────────────
                    # Female: Sample-level Cochran (QI 8b)
                    # ─────────────────────────────────────────────────────────────
                    if preproc_flags.get("samp_cochran", True):
                        female_df_sc = female_df.copy()
                        female_df_sc = _apply_cochran_selections(female_df_sc, st.session_state.get("rep_outlier_selections_female", {}))
                        female_df_sc = _apply_cochran_selections(female_df_sc, st.session_state.get("rep_bartlett_selections_female", {}))

                        female_samp_outliers = _detect_sample_cochran_outliers(female_df_sc)

                        if female_samp_outliers:
                            if st.session_state.get("outlier_mode", "auto") == "auto":
                                for out in female_samp_outliers:
                                    st.session_state["samp_cochran_selections_female"][(out["Subject"], out["Sample"])] = "remove"
                            else:
                                st.markdown("#### 🔬 Female: Sample-level Outliers (Cochran QI 8b)")
                                st.info(f"**Iterative Process:** Found **{len(female_samp_outliers)}** outlier sample(s).")

                                for i, outlier in enumerate(female_samp_outliers):
                                    subj = outlier["Subject"]
                                    samp = outlier["Sample"]
                                    variance = outlier["variance"]
                                    key = (subj, samp)
                                    step_num = i + 1

                                    with st.container():
                                        st.markdown(
                                            f"<div style='background: linear-gradient(90deg, #d946ef22 0%, transparent 100%); "
                                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #d946ef; margin-bottom: 0.5rem;'>"
                                            f"<b>Step {step_num}</b> — Subject {subj}, Sample {samp}</div>",
                                            unsafe_allow_html=True
                                        )
                                        col1, col2 = st.columns([1, 2])
                                        with col1:
                                            st.caption(f"**Variance:** s² = {variance:.4f}")
                                            st.caption(f"**G:** {outlier['G']:.3f} (Crit: {outlier['G_crit']:.3f})")
                                        with col2:
                                            options = ["🗑️ Remove Sample", "⏩ Ignore (Stop)"]
                                            stored_idx = 0
                                            if key in st.session_state.get("samp_cochran_selections_female", {}):
                                                if st.session_state["samp_cochran_selections_female"][key] == "ignore":
                                                    stored_idx = 1

                                            selected = st.radio("Action:", options, index=stored_idx,
                                                              key=f"female_samp_coch_{i}_{subj}_{samp}", horizontal=True, label_visibility="collapsed")

                                            action = "remove" if selected.startswith("🗑️") else "ignore"
                                            st.session_state["samp_cochran_selections_female"][key] = action

                                    if i < len(female_samp_outliers) - 1:
                                        st.markdown("<div style='text-align: center; color: #d946ef; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)
                                st.markdown("---")
                        else:
                            st.session_state["samp_cochran_selections_female"] = {}

                    # ─────────────────────────────────────────────────────────────
                    # Female: Within-Subject Bartlett (QI 10)
                    # ─────────────────────────────────────────────────────────────
                    if preproc_flags.get("wp_bartlett_exclusion", False):
                        female_df_wb = female_df.copy()
                        female_df_wb = _apply_cochran_selections(female_df_wb, st.session_state.get("rep_outlier_selections_female", {}))
                        female_df_wb = _apply_cochran_selections(female_df_wb, st.session_state.get("rep_bartlett_selections_female", {}))
                        for (subj, samp), action in st.session_state.get("samp_cochran_selections_female", {}).items():
                            if action == "remove":
                                female_df_wb = female_df_wb[~((female_df_wb["Subject"] == subj) & (female_df_wb["Sample"] == samp))]

                        female_wb_outliers = _detect_wp_bartlett_outliers(female_df_wb)

                        if female_wb_outliers:
                            if st.session_state.get("outlier_mode", "auto") == "auto":
                                for out in female_wb_outliers:
                                    st.session_state["wp_bartlett_selections_female"][out["Subject"]] = "remove"
                            else:
                                st.markdown("#### 🧬 Female: Within-Subject Variance Outliers (QI 10)")
                                st.info(f"**Iterative Process:** Found **{len(female_wb_outliers)}** subject(s) with heterogeneous variance.")

                                for i, outlier in enumerate(female_wb_outliers):
                                    subj = outlier["Subject"]
                                    variance = outlier["within_subject_variance"]
                                    sample_means = outlier["sample_means"]
                                    step_num = i + 1

                                    with st.container():
                                        st.markdown(
                                            f"<div style='background: linear-gradient(90deg, #ec489922 0%, transparent 100%); "
                                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #ec4899; margin-bottom: 0.5rem;'>"
                                            f"<b>Step {step_num}</b> — Subject {subj}</div>",
                                            unsafe_allow_html=True
                                        )
                                        col1, col2 = st.columns([1, 2])
                                        with col1:
                                            st.caption(f"**Within-subject s²:** {variance:.4f}")
                                            st.caption(f"**Bartlett p:** {outlier['p']:.4f}")
                                            means_str = " · ".join([f"S{sm['Sample']}={sm['Result']:.2f}" for sm in sample_means])
                                            st.code(means_str, language=None)

                                        with col2:
                                            options = [f"🗑️ Remove subject (proceed to Step {step_num + 1})", "⏩ Ignore (accept heterogeneity)"]
                                            stored_idx = 0
                                            if subj in st.session_state.get("wp_bartlett_selections_female", {}):
                                                stored_action = st.session_state["wp_bartlett_selections_female"][subj]
                                                if stored_action == "remove": stored_idx = 0
                                                elif stored_action == "ignore": stored_idx = 1

                                            selected = st.radio("Action:", options, index=stored_idx,
                                                              key=f"female_wp_bart_{i}_{subj}", horizontal=True, label_visibility="collapsed")

                                            action = "remove" if selected.startswith("🗑️") else "ignore"
                                            st.session_state["wp_bartlett_selections_female"][subj] = action

                                    if i < len(female_wb_outliers) - 1:
                                        st.markdown("<div style='text-align: center; color: #ec4899; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)

                                st.markdown(
                                    "<div style='background: #22c55e22; padding: 0.5rem; border-left: 3px solid #22c55e; margin-top: 0.5rem;'>"
                                    "✅ <b>Homogeneity Target</b> — After processing all steps, within-subject variances should be homogeneous.</div>",
                                    unsafe_allow_html=True
                                )
                                st.markdown("---")
                        else:
                            st.session_state["wp_bartlett_selections_female"] = {}

                    # ─────────────────────────────────────────────────────────────
                    # Female: Reed (QI 8c)
                    # ─────────────────────────────────────────────────────────────
                    if preproc_flags.get("reed", True):
                        female_df_reed = female_df.copy()
                        female_df_reed = _apply_cochran_selections(female_df_reed, st.session_state.get("rep_outlier_selections_female", {}))
                        female_df_reed = _apply_cochran_selections(female_df_reed, st.session_state.get("rep_bartlett_selections_female", {}))
                        for (subj, samp), action in st.session_state.get("samp_cochran_selections_female", {}).items():
                            if action == "remove":
                                female_df_reed = female_df_reed[~((female_df_reed["Subject"] == subj) & (female_df_reed["Sample"] == samp))]
                        for subj, action in st.session_state.get("wp_bartlett_selections_female", {}).items():
                            if action == "remove":
                                female_df_reed = female_df_reed[female_df_reed["Subject"] != subj]

                        female_reed_outliers = _detect_reed_outliers(female_df_reed)

                        if female_reed_outliers:
                            if st.session_state.get("outlier_mode", "auto") == "auto":
                                for out in female_reed_outliers:
                                    st.session_state["reed_selections_female"][out["Subject"]] = "remove"
                            else:
                                st.markdown("#### 🔬 Female: Reed Outliers (QI 8c)")
                                st.info(f"**Between-subject check:** Found **{len(female_reed_outliers)}** outlier subject(s).")

                                for i, outlier in enumerate(female_reed_outliers):
                                    subj = outlier["Subject"]
                                    mean_val = outlier["mean"]
                                    reason = outlier["reason"]
                                    step_num = i + 1

                                    with st.container():
                                        st.markdown(
                                            f"<div style='background: linear-gradient(90deg, #ef444422 0%, transparent 100%); "
                                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #ef4444; margin-bottom: 0.5rem;'>"
                                            f"<b>Step {step_num}</b> — Subject {subj}</div>",
                                            unsafe_allow_html=True
                                        )
                                        col1, col2 = st.columns([1, 2])
                                        with col1:
                                            st.caption(f"**Mean:** {mean_val:.4f}")
                                            st.caption(f"**Reason:** {reason}")
                                        with col2:
                                            options = ["🗑️ Remove Subject", "⏩ Ignore (Stop)"]
                                            stored_idx = 0
                                            if subj in st.session_state.get("reed_selections_female", {}):
                                                if st.session_state["reed_selections_female"][subj] == "ignore":
                                                    stored_idx = 1

                                            selected = st.radio("Action:", options, index=stored_idx,
                                                              key=f"female_reed_{i}_{subj}", horizontal=True, label_visibility="collapsed")

                                            action = "remove" if selected.startswith("🗑️") else "ignore"
                                            st.session_state["reed_selections_female"][subj] = action

                                    if i < len(female_reed_outliers) - 1:
                                        st.markdown("<div style='text-align: center; color: #ef4444; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)
                                st.markdown("---")
                        else:
                            if st.session_state.get("outlier_mode", "auto") != "auto":
                                st.markdown("#### 🔬 Female: Reed Outliers (QI 8c)")
                                st.success("✅ No outliers detected (or masked by prior exclusion steps).")
                                st.markdown("---")
                            st.session_state["reed_selections_female"] = {}

                    # ─────────────────────────────────────────────────────────────
                    # Female: Steady-State Drift (QI 7)
                    # ─────────────────────────────────────────────────────────────
                    if preproc_flags.get("drift", True):
                        female_df_drift = female_df.copy()
                        female_df_drift = _apply_cochran_selections(female_df_drift, st.session_state.get("rep_outlier_selections_female", {}))
                        female_df_drift = _apply_cochran_selections(female_df_drift, st.session_state.get("rep_bartlett_selections_female", {}))
                        for (subj, samp), action in st.session_state.get("samp_cochran_selections_female", {}).items():
                            if action == "remove":
                                female_df_drift = female_df_drift[~((female_df_drift["Subject"] == subj) & (female_df_drift["Sample"] == samp))]
                        for subj, action in st.session_state.get("wp_bartlett_selections_female", {}).items():
                            if action == "remove":
                                female_df_drift = female_df_drift[female_df_drift["Subject"] != subj]
                        for subj, action in st.session_state.get("reed_selections_female", {}).items():
                            if action == "remove":
                                female_df_drift = female_df_drift[female_df_drift["Subject"] != subj]

                        female_drift_outliers = _detect_drift_outliers(female_df_drift)

                        if female_drift_outliers:
                            if st.session_state.get("outlier_mode", "auto") == "auto":
                                for out in female_drift_outliers:
                                    st.session_state["drift_selections_female"][out["Subject"]] = "remove"
                            else:
                                st.markdown("#### 📉 Female: Steady-State Drift (QI 7)")
                                st.info(f"**Trend check:** Found **{len(female_drift_outliers)}** subject(s) with significant drift.")

                                for i, outlier in enumerate(female_drift_outliers):
                                    subj = outlier["Subject"]
                                    slope = outlier["slope"]
                                    p_val = outlier["p"]
                                    step_num = i + 1

                                    with st.container():
                                        st.markdown(
                                            f"<div style='background: linear-gradient(90deg, #d946ef22 0%, transparent 100%); "
                                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #d946ef; margin-bottom: 0.5rem;'>"
                                            f"<b>Step {step_num}</b> — Subject {subj}</div>",
                                            unsafe_allow_html=True
                                        )
                                        col1, col2 = st.columns([1, 2])
                                        with col1:
                                            st.caption(f"**Slope:** {slope:.4f}")
                                            st.caption(f"**p-value:** {p_val:.4f}")
                                        with col2:
                                            options = ["🗑️ Remove Subject", "⏩ Ignore (Accept Trend)"]
                                            stored_idx = 0
                                            if subj in st.session_state.get("drift_selections_female", {}):
                                                if st.session_state["drift_selections_female"][subj] == "ignore":
                                                    stored_idx = 1

                                            selected = st.radio("Action:", options, index=stored_idx,
                                                              key=f"female_drift_{i}_{subj}", horizontal=True, label_visibility="collapsed")

                                            action = "remove" if selected.startswith("🗑️") else "ignore"
                                            st.session_state["drift_selections_female"][subj] = action

                                    if i < len(female_drift_outliers) - 1:
                                        st.markdown("<div style='text-align: center; color: #d946ef; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)
                                st.markdown("---")
                        else:
                            st.session_state["drift_selections_female"] = {}

                # ═══════════════════════════════════════════════════════════════
                # MALE - All Outlier Steps
                # ═══════════════════════════════════════════════════════════════
                if not male_df.empty:
                    if st.session_state.get("outlier_mode", "auto") != "auto":
                        st.markdown("---")
                        st.markdown("## 👨 Male Outlier Analysis")
                        st.markdown("---")

                    # ─────────────────────────────────────────────────────────────
                    # Male: Replicate-level Cochran (QI 8a)
                    # ─────────────────────────────────────────────────────────────
                    male_outliers = _detect_replicate_outliers(male_df)
                    if male_outliers:
                        if st.session_state.get("outlier_mode", "auto") == "auto":
                            for outlier in male_outliers:
                                key = (outlier["Subject"], outlier["Sample"])
                                st.session_state["rep_outlier_selections_male"][key] = "remove_all"
                        else:
                            st.markdown("#### 🔍 Male: Replicate-level Outliers (Cochran QI 8a)")
                            st.caption(f"Found **{len(male_outliers)}** sample(s) with high variance.")

                            for i, outlier in enumerate(male_outliers):
                                subj, samp = outlier["Subject"], outlier["Sample"]
                                variance, reps = outlier["variance"], outlier["replicates"]
                                key = (subj, samp)

                                with st.container():
                                    col1, col2 = st.columns([1, 2])

                                    with col1:
                                        st.markdown(f"**Subject {subj}, Sample {samp}**")
                                        st.caption(f"s² = {variance:.4f} | G = {outlier['G']:.3f}")
                                        rep_str = " · ".join([f"R{r['Replicate']}={r['Result']:.2f}" for r in reps])
                                        st.code(rep_str, language=None)

                                    with col2:
                                        options = ["🗑️ Remove sample"]
                                        for r in reps:
                                            options.append(f"✅ Keep R{r['Replicate']}")
                                        options.append("⏩ Ignore")

                                        stored_idx = 0
                                        if key in st.session_state.get("rep_outlier_selections_male", {}):
                                            stored_action = st.session_state["rep_outlier_selections_male"][key]
                                            if stored_action == "remove_all": stored_idx = 0
                                            elif stored_action == "ignore": stored_idx = len(options) - 1
                                            elif stored_action.startswith("keep_rep_"):
                                                for j, r in enumerate(reps):
                                                    if str(r["Replicate"]) == stored_action.replace("keep_rep_", ""):
                                                        stored_idx = j + 1; break

                                        selected = st.radio("Action:", options, index=stored_idx,
                                            key=f"male_outlier_{i}_{subj}_{samp}", horizontal=True, label_visibility="collapsed")

                                        if selected.startswith("🗑️"): action = "remove_all"
                                        elif selected.startswith("⏩"): action = "ignore"
                                        else:
                                            for r in reps:
                                                if f"R{r['Replicate']}" in selected:
                                                    action = f"keep_rep_{r['Replicate']}"; break
                                        st.session_state["rep_outlier_selections_male"][key] = action

                                if i < len(male_outliers) - 1:
                                    st.markdown("<hr style='margin: 0.5rem 0; opacity: 0.3;'>", unsafe_allow_html=True)
                            st.markdown("---")
                    else:
                        st.session_state["rep_outlier_selections_male"] = {}

                    # ─────────────────────────────────────────────────────────────
                    # Male: Replicate-level Bartlett (QI 8a)
                    # ─────────────────────────────────────────────────────────────
                    if preproc_flags.get("rep_bartlett_exclusion", False):
                        male_df_bartlett = _apply_cochran_selections(male_df.copy(), st.session_state.get("rep_outlier_selections_male", {}))
                        male_bartlett_outliers = _detect_bartlett_outliers(male_df_bartlett)

                        if male_bartlett_outliers:
                            if st.session_state.get("outlier_mode", "auto") == "auto":
                                for outlier in male_bartlett_outliers:
                                    key = (outlier["Subject"], outlier["Sample"])
                                    st.session_state["rep_bartlett_selections_male"][key] = "remove_all"
                            else:
                                st.markdown("#### 🔬 Male: Bartlett Iterative Outliers (QI 8a)")
                                st.info(f"**After Cochran cleanup:** Found **{len(male_bartlett_outliers)}** sample(s) for iterative Bartlett exclusion.")

                                for i, outlier in enumerate(male_bartlett_outliers):
                                    subj, samp = outlier["Subject"], outlier["Sample"]
                                    variance, reps = outlier["variance"], outlier["replicates"]
                                    key = (subj, samp)
                                    step_num = i + 1

                                    with st.container():
                                        st.markdown(
                                            f"<div style='background: linear-gradient(90deg, #3b82f622 0%, transparent 100%); "
                                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #3b82f6; margin-bottom: 0.5rem;'>"
                                            f"<b>Step {step_num}</b> — Subject {subj}, Sample {samp}</div>",
                                            unsafe_allow_html=True
                                        )
                                        col1, col2 = st.columns([1, 2])

                                        with col1:
                                            st.caption(f"**Variance:** s² = {variance:.4f}")
                                            rep_str = " · ".join([f"R{r['Replicate']}={r['Result']:.2f}" for r in reps])
                                            st.code(rep_str, language=None)

                                        with col2:
                                            options = [f"🗑️ Remove (proceed to Step {step_num + 1})"]
                                            for r in reps:
                                                options.append(f"✅ Keep only R{r['Replicate']}")
                                            options.append("⏩ Ignore (accept heterogeneity)")

                                            stored_idx = 0
                                            if key in st.session_state.get("rep_bartlett_selections_male", {}):
                                                stored_action = st.session_state["rep_bartlett_selections_male"][key]
                                                if stored_action == "remove_all": stored_idx = 0
                                                elif stored_action == "ignore": stored_idx = len(options) - 1
                                                elif stored_action.startswith("keep_rep_"):
                                                    for j, r in enumerate(reps):
                                                        if str(r["Replicate"]) == stored_action.replace("keep_rep_", ""):
                                                            stored_idx = j + 1; break

                                            selected = st.radio("Action:", options, index=stored_idx,
                                                key=f"male_bartlett_{i}_{subj}_{samp}", horizontal=True, label_visibility="collapsed")

                                            if selected.startswith("🗑️"): action = "remove_all"
                                            elif selected.startswith("⏩"): action = "ignore"
                                            else:
                                                for r in reps:
                                                    if f"R{r['Replicate']}" in selected:
                                                        action = f"keep_rep_{r['Replicate']}"; break
                                            st.session_state["rep_bartlett_selections_male"][key] = action

                                    if i < len(male_bartlett_outliers) - 1:
                                        st.markdown("<div style='text-align: center; color: #3b82f6; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)
                                st.markdown(
                                    "<div style='background: #22c55e22; padding: 0.5rem; border-left: 3px solid #22c55e; margin-top: 0.5rem;'>"
                                    "✅ <b>Homogeneity Target</b> — After processing all steps above, variances should be homogeneous.</div>",
                                    unsafe_allow_html=True
                                )
                                st.markdown("---")
                        else:
                            st.session_state["rep_bartlett_selections_male"] = {}

                    # ─────────────────────────────────────────────────────────────
                    # Male: Sample-level Cochran (QI 8b)
                    # ─────────────────────────────────────────────────────────────
                    if preproc_flags.get("samp_cochran", True):
                        male_df_sc = male_df.copy()
                        male_df_sc = _apply_cochran_selections(male_df_sc, st.session_state.get("rep_outlier_selections_male", {}))
                        male_df_sc = _apply_cochran_selections(male_df_sc, st.session_state.get("rep_bartlett_selections_male", {}))

                        male_samp_outliers = _detect_sample_cochran_outliers(male_df_sc)

                        if male_samp_outliers:
                            if st.session_state.get("outlier_mode", "auto") == "auto":
                                for out in male_samp_outliers:
                                    st.session_state["samp_cochran_selections_male"][(out["Subject"], out["Sample"])] = "remove"
                            else:
                                st.markdown("#### 🔬 Male: Sample-level Outliers (Cochran QI 8b)")
                                st.info(f"**Iterative Process:** Found **{len(male_samp_outliers)}** outlier sample(s).")

                                for i, outlier in enumerate(male_samp_outliers):
                                    subj = outlier["Subject"]
                                    samp = outlier["Sample"]
                                    variance = outlier["variance"]
                                    key = (subj, samp)
                                    step_num = i + 1

                                    with st.container():
                                        st.markdown(
                                            f"<div style='background: linear-gradient(90deg, #8b5cf622 0%, transparent 100%); "
                                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #8b5cf6; margin-bottom: 0.5rem;'>"
                                            f"<b>Step {step_num}</b> — Subject {subj}, Sample {samp}</div>",
                                            unsafe_allow_html=True
                                        )
                                        col1, col2 = st.columns([1, 2])
                                        with col1:
                                            st.caption(f"**Variance:** s² = {variance:.4f}")
                                            st.caption(f"**G:** {outlier['G']:.3f} (Crit: {outlier['G_crit']:.3f})")
                                        with col2:
                                            options = ["🗑️ Remove Sample", "⏩ Ignore (Stop)"]
                                            stored_idx = 0
                                            if key in st.session_state.get("samp_cochran_selections_male", {}):
                                                if st.session_state["samp_cochran_selections_male"][key] == "ignore":
                                                    stored_idx = 1

                                            selected = st.radio("Action:", options, index=stored_idx,
                                                              key=f"male_samp_coch_{i}_{subj}_{samp}", horizontal=True, label_visibility="collapsed")

                                            action = "remove" if selected.startswith("🗑️") else "ignore"
                                            st.session_state["samp_cochran_selections_male"][key] = action

                                    if i < len(male_samp_outliers) - 1:
                                        st.markdown("<div style='text-align: center; color: #8b5cf6; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)
                                st.markdown("---")
                        else:
                            st.session_state["samp_cochran_selections_male"] = {}

                    # ─────────────────────────────────────────────────────────────
                    # Male: Within-Subject Bartlett (QI 10)
                    # ─────────────────────────────────────────────────────────────
                    if preproc_flags.get("wp_bartlett_exclusion", False):
                        male_df_wb = male_df.copy()
                        male_df_wb = _apply_cochran_selections(male_df_wb, st.session_state.get("rep_outlier_selections_male", {}))
                        male_df_wb = _apply_cochran_selections(male_df_wb, st.session_state.get("rep_bartlett_selections_male", {}))
                        for (subj, samp), action in st.session_state.get("samp_cochran_selections_male", {}).items():
                            if action == "remove":
                                male_df_wb = male_df_wb[~((male_df_wb["Subject"] == subj) & (male_df_wb["Sample"] == samp))]

                        male_wb_outliers = _detect_wp_bartlett_outliers(male_df_wb)

                        if male_wb_outliers:
                            if st.session_state.get("outlier_mode", "auto") == "auto":
                                for out in male_wb_outliers:
                                    st.session_state["wp_bartlett_selections_male"][out["Subject"]] = "remove"
                            else:
                                st.markdown("#### 🧬 Male: Within-Subject Variance Outliers (QI 10)")
                                st.info(f"**Iterative Process:** Found **{len(male_wb_outliers)}** subject(s) with heterogeneous variance.")

                                for i, outlier in enumerate(male_wb_outliers):
                                    subj = outlier["Subject"]
                                    variance = outlier["within_subject_variance"]
                                    sample_means = outlier["sample_means"]
                                    step_num = i + 1

                                    with st.container():
                                        st.markdown(
                                            f"<div style='background: linear-gradient(90deg, #f59e0b22 0%, transparent 100%); "
                                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #f59e0b; margin-bottom: 0.5rem;'>"
                                            f"<b>Step {step_num}</b> — Subject {subj}</div>",
                                            unsafe_allow_html=True
                                        )
                                        col1, col2 = st.columns([1, 2])
                                        with col1:
                                            st.caption(f"**Within-subject s²:** {variance:.4f}")
                                            st.caption(f"**Bartlett p:** {outlier['p']:.4f}")
                                            means_str = " · ".join([f"S{sm['Sample']}={sm['Result']:.2f}" for sm in sample_means])
                                            st.code(means_str, language=None)

                                        with col2:
                                            options = [f"🗑️ Remove subject (proceed to Step {step_num + 1})", "⏩ Ignore (accept heterogeneity)"]
                                            stored_idx = 0
                                            if subj in st.session_state.get("wp_bartlett_selections_male", {}):
                                                stored_action = st.session_state["wp_bartlett_selections_male"][subj]
                                                if stored_action == "remove": stored_idx = 0
                                                elif stored_action == "ignore": stored_idx = 1

                                            selected = st.radio("Action:", options, index=stored_idx,
                                                              key=f"male_wp_bart_{i}_{subj}", horizontal=True, label_visibility="collapsed")

                                            action = "remove" if selected.startswith("🗑️") else "ignore"
                                            st.session_state["wp_bartlett_selections_male"][subj] = action

                                    if i < len(male_wb_outliers) - 1:
                                        st.markdown("<div style='text-align: center; color: #f59e0b; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)

                                st.markdown(
                                    "<div style='background: #22c55e22; padding: 0.5rem; border-left: 3px solid #22c55e; margin-top: 0.5rem;'>"
                                    "✅ <b>Homogeneity Target</b> — After processing all steps, within-subject variances should be homogeneous.</div>",
                                    unsafe_allow_html=True
                                )
                                st.markdown("---")
                        else:
                            st.session_state["wp_bartlett_selections_male"] = {}

                    # ─────────────────────────────────────────────────────────────
                    # Male: Reed (QI 8c)
                    # ─────────────────────────────────────────────────────────────
                    if preproc_flags.get("reed", True):
                        male_df_reed = male_df.copy()
                        male_df_reed = _apply_cochran_selections(male_df_reed, st.session_state.get("rep_outlier_selections_male", {}))
                        male_df_reed = _apply_cochran_selections(male_df_reed, st.session_state.get("rep_bartlett_selections_male", {}))
                        for (subj, samp), action in st.session_state.get("samp_cochran_selections_male", {}).items():
                            if action == "remove":
                                male_df_reed = male_df_reed[~((male_df_reed["Subject"] == subj) & (male_df_reed["Sample"] == samp))]
                        for subj, action in st.session_state.get("wp_bartlett_selections_male", {}).items():
                            if action == "remove":
                                male_df_reed = male_df_reed[male_df_reed["Subject"] != subj]

                        male_reed_outliers = _detect_reed_outliers(male_df_reed)

                        if male_reed_outliers:
                            if st.session_state.get("outlier_mode", "auto") == "auto":
                                for out in male_reed_outliers:
                                    st.session_state["reed_selections_male"][out["Subject"]] = "remove"
                            else:
                                st.markdown("#### 🔬 Male: Reed Outliers (QI 8c)")
                                st.info(f"**Between-subject check:** Found **{len(male_reed_outliers)}** outlier subject(s).")

                                for i, outlier in enumerate(male_reed_outliers):
                                    subj = outlier["Subject"]
                                    mean_val = outlier["mean"]
                                    reason = outlier["reason"]
                                    step_num = i + 1

                                    with st.container():
                                        st.markdown(
                                            f"<div style='background: linear-gradient(90deg, #ef444422 0%, transparent 100%); "
                                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #ef4444; margin-bottom: 0.5rem;'>"
                                            f"<b>Step {step_num}</b> — Subject {subj}</div>",
                                            unsafe_allow_html=True
                                        )
                                        col1, col2 = st.columns([1, 2])
                                        with col1:
                                            st.caption(f"**Mean:** {mean_val:.4f}")
                                            st.caption(f"**Reason:** {reason}")
                                        with col2:
                                            options = ["🗑️ Remove Subject", "⏩ Ignore (Stop)"]
                                            stored_idx = 0
                                            if subj in st.session_state.get("reed_selections_male", {}):
                                                if st.session_state["reed_selections_male"][subj] == "ignore":
                                                    stored_idx = 1

                                            selected = st.radio("Action:", options, index=stored_idx,
                                                              key=f"male_reed_{i}_{subj}", horizontal=True, label_visibility="collapsed")

                                            action = "remove" if selected.startswith("🗑️") else "ignore"
                                            st.session_state["reed_selections_male"][subj] = action

                                    if i < len(male_reed_outliers) - 1:
                                        st.markdown("<div style='text-align: center; color: #ef4444; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)
                                st.markdown("---")
                        else:
                            if st.session_state.get("outlier_mode", "auto") != "auto":
                                st.markdown("#### 🔬 Male: Reed Outliers (QI 8c)")
                                st.success("✅ No outliers detected (or masked by prior exclusion steps).")
                                st.markdown("---")
                            st.session_state["reed_selections_male"] = {}

                    # ─────────────────────────────────────────────────────────────
                    # Male: Steady-State Drift (QI 7)
                    # ─────────────────────────────────────────────────────────────
                    if preproc_flags.get("drift", True):
                        male_df_drift = male_df.copy()
                        male_df_drift = _apply_cochran_selections(male_df_drift, st.session_state.get("rep_outlier_selections_male", {}))
                        male_df_drift = _apply_cochran_selections(male_df_drift, st.session_state.get("rep_bartlett_selections_male", {}))
                        for (subj, samp), action in st.session_state.get("samp_cochran_selections_male", {}).items():
                            if action == "remove":
                                male_df_drift = male_df_drift[~((male_df_drift["Subject"] == subj) & (male_df_drift["Sample"] == samp))]
                        for subj, action in st.session_state.get("wp_bartlett_selections_male", {}).items():
                            if action == "remove":
                                male_df_drift = male_df_drift[male_df_drift["Subject"] != subj]
                        for subj, action in st.session_state.get("reed_selections_male", {}).items():
                            if action == "remove":
                                male_df_drift = male_df_drift[male_df_drift["Subject"] != subj]

                        male_drift_outliers = _detect_drift_outliers(male_df_drift)

                        if male_drift_outliers:
                            if st.session_state.get("outlier_mode", "auto") == "auto":
                                for out in male_drift_outliers:
                                    st.session_state["drift_selections_male"][out["Subject"]] = "remove"
                            else:
                                st.markdown("#### 📉 Male: Steady-State Drift (QI 7)")
                                st.info(f"**Trend check:** Found **{len(male_drift_outliers)}** subject(s) with significant drift.")

                                for i, outlier in enumerate(male_drift_outliers):
                                    subj = outlier["Subject"]
                                    slope = outlier["slope"]
                                    p_val = outlier["p"]
                                    step_num = i + 1

                                    with st.container():
                                        st.markdown(
                                            f"<div style='background: linear-gradient(90deg, #8b5cf622 0%, transparent 100%); "
                                            f"padding: 0.3rem 0.5rem; border-left: 3px solid #8b5cf6; margin-bottom: 0.5rem;'>"
                                            f"<b>Step {step_num}</b> — Subject {subj}</div>",
                                            unsafe_allow_html=True
                                        )
                                        col1, col2 = st.columns([1, 2])
                                        with col1:
                                            st.caption(f"**Slope:** {slope:.4f}")
                                            st.caption(f"**p-value:** {p_val:.4f}")
                                        with col2:
                                            options = ["🗑️ Remove Subject", "⏩ Ignore (Accept Trend)"]
                                            stored_idx = 0
                                            if subj in st.session_state.get("drift_selections_male", {}):
                                                if st.session_state["drift_selections_male"][subj] == "ignore":
                                                    stored_idx = 1

                                            selected = st.radio("Action:", options, index=stored_idx,
                                                              key=f"male_drift_{i}_{subj}", horizontal=True, label_visibility="collapsed")

                                            action = "remove" if selected.startswith("🗑️") else "ignore"
                                            st.session_state["drift_selections_male"][subj] = action

                                    if i < len(male_drift_outliers) - 1:
                                        st.markdown("<div style='text-align: center; color: #8b5cf6; font-size: 1.5rem; margin: 0.3rem 0;'>↓</div>", unsafe_allow_html=True)
                                st.markdown("---")
                        else:
                            st.session_state["drift_selections_male"] = {}
    else:
        # Clear gender-specific selections if not applicable
        st.session_state["rep_outlier_selections_male"] = {}
        st.session_state["rep_outlier_selections_female"] = {}
        st.session_state["rep_bartlett_selections_male"] = {}
        st.session_state["rep_bartlett_selections_female"] = {}
        st.session_state["samp_cochran_selections_male"] = {}
        st.session_state["samp_cochran_selections_female"] = {}
        st.session_state["wp_bartlett_selections_male"] = {}
        st.session_state["wp_bartlett_selections_female"] = {}
        st.session_state["reed_selections_male"] = {}
        st.session_state["reed_selections_female"] = {}
        st.session_state["drift_selections_male"] = {}
        st.session_state["drift_selections_female"] = {}

    calc_disabled = (not st.session_state.get("mapping_ok", False)) or (
        st.session_state.get("gender_based", False) and (not st.session_state.get("gender_settings_ok", True))
    )

    if calc_disabled:
        st.info("Save a valid column mapping above to enable **Calculate**.")

    # CHANGE: add disabled=calc_disabled
    if st.button("Calculate", type="primary", disabled=calc_disabled):
        try:
            df_for_calc = st.session_state.get("mapped_df")
            # Keep the extra guard anyway
            if df_for_calc is None:
                st.warning("Please open “Map / confirm columns” and save a valid mapping first.")
                st.stop()

            try:
                res = calculate_bv(
                    df_for_calc,
                    use_cv_anova=st.session_state.get("use_cv_anova", False),
                    enforce_balance=st.session_state.get("enforce_balance", True),
                    rep_outlier_selections=st.session_state.get("rep_outlier_selections", {}),
                    rep_bartlett_selections=st.session_state.get("rep_bartlett_selections", {}),
                    samp_cochran_selections=st.session_state.get("samp_cochran_selections", {}),
                    wp_bartlett_selections=st.session_state.get("wp_bartlett_selections", {}),
                    reed_selections=st.session_state.get("reed_selections", {}),
                    drift_selections=st.session_state.get("drift_selections", {}),
                    mapping_warnings=st.session_state.get("mapping_warnings", []),
                )
                # ─────────────────────────────────────────────────────────────
                # --- Gender-split calculations (run the same pipeline per gender) ---
                res_male = res_female = None
                male_df_final = female_df_final = None

                if gender_based:
                    df_gender = df_for_calc.copy()

                    mode = st.session_state.get("gender_mode", "Single gender column")
                    if mode == "Single gender column":
                        gcol = st.session_state["gender_col"]
                        mv = st.session_state["male_value"]
                        fv = st.session_state["female_value"]
                        df_gender["Gender"] = _derive_gender_single_column(user_df, gcol, mv, fv)
                    else:
                        mcol = st.session_state["male_col"]
                        fcol = st.session_state["female_col"]
                        df_gender["Gender"] = _derive_gender_two_indicators(user_df, mcol, fcol)

                    # keep only clear Male/Female rows
                    df_gender = df_gender[df_gender["Gender"].isin(["Male", "Female"])].copy()

                    # enforce one gender per subject
                    df_gender, g_log = _enforce_subject_level_gender(df_gender)
                    for ln in g_log:
                        st.info(ln)

                    n_m = int((df_gender["Gender"] == "Male").sum())
                    n_f = int((df_gender["Gender"] == "Female").sum())
                    st.divider()
                    st.write(f"**Gender split rows available:** Male = {n_m}, Female = {n_f}")

                    if (n_m == 0) or (n_f == 0):
                        st.error("Gender-based calculations enabled, but one gender group has 0 usable rows after mapping.")
                    else:
                        use_cv_anova_flag = st.session_state.get("use_cv_anova", False)
                        enforce_balance_flag = st.session_state.get("enforce_balance", True)

                        # (optional but recommended) drop Gender column before BV calc
                        male_input = df_gender[df_gender["Gender"] == "Male"].drop(columns=["Gender"])
                        female_input = df_gender[df_gender["Gender"] == "Female"].drop(columns=["Gender"])

                        # Calculate gender-specific mapping warnings for dropped rows
                        dropped_indices = st.session_state.get("dropped_row_indices", [])
                        male_warnings = []
                        female_warnings = []
                        if dropped_indices:
                            # Determine gender for each dropped row from original user_df
                            mode = st.session_state.get("gender_mode", "Single gender column")
                            n_male_dropped = 0
                            n_female_dropped = 0
                            n_unknown_dropped = 0

                            for idx in dropped_indices:
                                if idx in user_df.index:
                                    if mode == "Single gender column":
                                        gcol = st.session_state.get("gender_col")
                                        mv = st.session_state.get("male_value")
                                        fv = st.session_state.get("female_value")
                                        if gcol and gcol in user_df.columns:
                                            val = user_df.loc[idx, gcol]
                                            if str(val) == str(mv):
                                                n_male_dropped += 1
                                            elif str(val) == str(fv):
                                                n_female_dropped += 1
                                            else:
                                                n_unknown_dropped += 1
                                    else:
                                        mcol = st.session_state.get("male_col")
                                        fcol = st.session_state.get("female_col")
                                        if mcol and fcol and mcol in user_df.columns and fcol in user_df.columns:
                                            m_val = user_df.loc[idx, mcol]
                                            f_val = user_df.loc[idx, fcol]
                                            is_m = pd.notna(m_val) and str(m_val).strip() not in ("", "0", "nan")
                                            is_f = pd.notna(f_val) and str(f_val).strip() not in ("", "0", "nan")
                                            if is_m and not is_f:
                                                n_male_dropped += 1
                                            elif is_f and not is_m:
                                                n_female_dropped += 1
                                            else:
                                                n_unknown_dropped += 1

                            if n_male_dropped > 0:
                                male_warnings.append(f"Dropped {n_male_dropped} row(s) with missing/empty Result values.")
                            if n_female_dropped > 0:
                                female_warnings.append(f"Dropped {n_female_dropped} row(s) with missing/empty Result values.")

                        try:
                            res_male = calculate_bv(
                                male_input,
                                use_cv_anova=use_cv_anova_flag,
                                enforce_balance=enforce_balance_flag,
                                rep_outlier_selections=st.session_state.get("rep_outlier_selections_male", {}),
                                rep_bartlett_selections=st.session_state.get("rep_bartlett_selections_male", {}),
                                samp_cochran_selections=st.session_state.get("samp_cochran_selections_male", {}),
                                wp_bartlett_selections=st.session_state.get("wp_bartlett_selections_male", {}),
                                reed_selections=st.session_state.get("reed_selections_male", {}),
                                drift_selections=st.session_state.get("drift_selections_male", {}),
                                mapping_warnings=male_warnings,
                            )
                            male_df_final = res_male.clean_df
                        except Exception as e:
                            st.error(f"Male calculation failed: {e}")

                        try:
                            res_female = calculate_bv(
                                female_input,
                                use_cv_anova=use_cv_anova_flag,
                                enforce_balance=enforce_balance_flag,
                                rep_outlier_selections=st.session_state.get("rep_outlier_selections_female", {}),
                                rep_bartlett_selections=st.session_state.get("rep_bartlett_selections_female", {}),
                                samp_cochran_selections=st.session_state.get("samp_cochran_selections_female", {}),
                                wp_bartlett_selections=st.session_state.get("wp_bartlett_selections_female", {}),
                                reed_selections=st.session_state.get("reed_selections_female", {}),
                                drift_selections=st.session_state.get("drift_selections_female", {}),
                                mapping_warnings=female_warnings,
                            )
                            female_df_final = res_female.clean_df
                        except Exception as e:
                            st.error(f"Female calculation failed: {e}")

                # ─────────────────────────────────────────────────────────────
                # Group registry for gender-based rendering
                # ─────────────────────────────────────────────────────────────
                groups: dict[str, tuple[BVResult, pd.DataFrame]] = {"Overall": (res, df_for_calc)}

                if st.session_state.get("gender_based", False) and (res_male is not None) and (res_female is not None):
                    groups["Male"] = (res_male, male_input)
                    groups["Female"] = (res_female, female_input)

                unit = st.session_state.get("result_unit", "").strip()
                use_cv_anova_flag = st.session_state.get("use_cv_anova", False)
                flags = st.session_state.get("preproc_flags", {})
                qi_manual = st.session_state.get("qi_manual")


                # ─────────────────────────────────────────────────────────────
                # Gender-based reporting tables (minimal, no rewrite of your UI cards)
                # ─────────────────────────────────────────────────────────────
                if st.session_state.get("gender_based", False) and (res_male is not None) and (res_female is not None):

                    use_cv_anova_flag = st.session_state.get("use_cv_anova", False)
                    unit = st.session_state.get("result_unit", "").strip()

                    def _one_row(label: str, r: BVResult) -> dict:
                        # --- Standard ANOVA CVI (always available) ---
                        cvi_std_val = float(r.cv_I)
                        cvi_std_ci  = tuple(r.ci_cv_I)

                        # --- CV-ANOVA CVI (optional) ---
                        if (r.cv_I_cv_anova is not None) and (r.ci_cv_I_cv_anova is not None):
                            cvi_cva_val = float(r.cv_I_cv_anova)
                            cvi_cva_ci  = tuple(r.ci_cv_I_cv_anova)
                            cvi_cva_txt = f"{cvi_cva_val:.2f} ({cvi_cva_ci[0]:.2f}–{cvi_cva_ci[1]:.2f})"
                        else:
                            cvi_cva_txt = "—"  # or "Not computed"

                        return {
                            "Gender": label,
                            "Mean concentration (95% CI)": f"{r.grand_mean:.2f}" + (f" {unit} " if unit else "") + f"({r.ci_mean[0]:.2f}–{r.ci_mean[1]:.2f})" + (f" {unit}" if unit else ""),
                            #"Mean CI (95%)": f"{r.ci_mean[0]:.2f}–{r.ci_mean[1]:.2f}" + (f" {unit}" if unit else ""),
                            "CVA % (95% CI)": f"{r.cv_A:.2f} ({r.ci_cv_A[0]:.2f}–{r.ci_cv_A[1]:.2f})",

                            # ✅ NEW: both CVIs as separate columns
                            "CVI (based on standard ANOVA) (95% CI)": f"{cvi_std_val:.2f} ({cvi_std_ci[0]:.2f}–{cvi_std_ci[1]:.2f})",
                            "CVI (based on CV-ANOVA) (95% CI)": cvi_cva_txt,

                            "CVG % (95% CI)": f"{r.cv_G:.2f} ({r.ci_cv_G[0]:.2f}–{r.ci_cv_G[1]:.2f})",
                            "RCV (95%)": f"−{r.rcv_95_down:.2f}% / +{r.rcv_95_up:.2f}%",
                        }



                    # ✅ NOW show recommendations AFTER the tables/plots
                    st.divider()
                    st.subheader("Gender-based reporting recommendations")

                    # Optional plots per gender (only if we successfully built final cleaned dfs)
                    tabs = st.tabs(["Male plot", "Female plot"])
                    with tabs[0]:
                        if male_df_final is not None and not male_df_final.empty:
                            st.plotly_chart(plot_subject_ranges(male_df_final), use_container_width=True)
                        else:
                            st.info("Male plot not available (cleaned dataset empty or failed).")
                    with tabs[1]:
                        if female_df_final is not None and not female_df_final.empty:
                            st.plotly_chart(plot_subject_ranges(female_df_final), use_container_width=True)
                        else:
                            st.info("Female plot not available (cleaned dataset empty or failed).")


                    #st.subheader("Gender-based BV results (side-by-side)")
                    gender_tbl = pd.DataFrame([
                        _one_row("Male", res_male),
                        _one_row("Female", res_female),
                    ])
                    # ✅ Rename columns to HTML subscript notation (uppercase A/I/G)
                    gender_tbl_disp = gender_tbl.rename(columns={
                        "CVA % (95% CI)": "CV<sub>A</sub> % (95% CI)",
                        "CVI (based on standard ANOVA) (95% CI)": "CV<sub>I</sub> (based on standard ANOVA) (95% CI)",
                        "CVI (based on CV-ANOVA) (95% CI)": "CV<sub>I</sub> (based on CV-ANOVA) (95% CI)",
                        "CVG % (95% CI)": "CV<sub>G</sub> % (95% CI)",
                    })

                    # ✅ Render as HTML so <sub> works in headers; hide index
                    gender_html = (
                        gender_tbl_disp.style
                            .set_table_styles([
                                {
                                    "selector": "th",
                                    "props": [
                                        ("background-color", "#d0e7f5"),
                                        ("color", "#1e3a5f"),
                                        ("font-weight", "700"),
                                        ("text-align", "left"),
                                    ],
                                },
                                {"selector": "td", "props": [("text-align", "left")]},
                            ])
                            .hide(axis="index")
                            .to_html(index=False, escape=False)   # escape=False allows <sub> to render
                    )

                    st.markdown(gender_html, unsafe_allow_html=True)

                    with st.expander("Gender-based reporting recommendations", expanded=True):
                        _render_gender_overlap_recommendations(res_male, res_female, use_cv_anova_flag)
                    st.divider()


                # ─────────────────────────────────────────────────────────────
                # ─────────────────────────────────────────────────────────────
                # Gender-aware render helpers (Key metrics / Summary / APS / BIVAC)
                # ─────────────────────────────────────────────────────────────

                def _fmt_ci(ci: tuple[float, float] | None, suffix: str = "") -> str:
                    if not ci:
                        return ""
                    lo, hi = float(ci[0]), float(ci[1])
                    return f"{lo:.2f}–{hi:.2f}{suffix}"

                def _fmt_count(x) -> str:
                    if x is None:
                        return ""
                    try:
                        xf = float(x)
                        return str(int(xf)) if xf.is_integer() else f"{xf:.2f}"
                    except Exception:
                        return str(x)

                def _pick_cvi_for_table(r: BVResult, use_cv_anova_flag: bool):
                    """Return (label, value, ci) for CVI row."""
                    if use_cv_anova_flag and (r.cv_I_cv_anova is not None) and (r.ci_cv_I_cv_anova is not None):
                        return "CVI (CV-ANOVA) (%)", float(r.cv_I_cv_anova), tuple(r.ci_cv_I_cv_anova)
                    return "CVI (standard ANOVA) (%)", float(r.cv_I), tuple(r.ci_cv_I)

                def _build_study_design_table(r: BVResult) -> pd.DataFrame:
                    # If unbalanced, S/R may be mean values (floats). Label accordingly.
                    try:
                        s_is_int = float(r.S).is_integer()
                        r_is_int = float(r.R).is_integer()
                    except Exception:
                        s_is_int = r_is_int = True

                    s_label = "Samples per subject (S)" if s_is_int else "Samples per subject (mean S̄)"
                    r_label = "Replicates per sample (R)" if r_is_int else "Replicates per sample (mean R̄)"

                    rows = [
                        {"Parameter": "Number of subjects (I)", "Value": _fmt_count(r.I)},
                        {"Parameter": s_label,                 "Value": _fmt_count(r.S)},
                        {"Parameter": r_label,                 "Value": _fmt_count(r.R)},
                    ]
                    return pd.DataFrame(rows)

                def _build_core_estimates_table(r: BVResult, unit: str, use_cv_anova_flag: bool) -> pd.DataFrame:
                    unit_txt = f" {unit}" if unit else ""

                    rows = [
                        {
                            "Parameter": f"Mean concentration {unit_txt}",
                            "Value": f"{r.grand_mean:.2f} {unit_txt}",
                            "CI (95%)": _fmt_ci(r.ci_mean, unit_txt),
                        },
                        {
                            "Parameter": "Analytical CV (%)",
                            "Value": f"{r.cv_A:.2f} %",
                            "CI (95%)": _fmt_ci(r.ci_cv_A, " %"),
                        },

                        # ✅ Always show standard ANOVA CVI
                        {
                            "Parameter": "Within-subject CV (based on standard ANOVA) (%)",
                            "Value": f"{float(r.cv_I):.2f} %",
                            "CI (95%)": _fmt_ci(tuple(r.ci_cv_I), " %"),
                        },
                    ]

                    # ✅ Also show CV-ANOVA CVI (if computed), otherwise show placeholder
                    if (r.cv_I_cv_anova is not None) and (r.ci_cv_I_cv_anova is not None):
                        rows.append({
                            "Parameter": "Within-subject CV (based on CV-ANOVA) (%)",
                            "Value": f"{float(r.cv_I_cv_anova):.2f} %",
                            "CI (95%)": _fmt_ci(tuple(r.ci_cv_I_cv_anova), " %"),
                        })
                    else:
                        rows.append({
                            "Parameter": "Within-subject CV (based on CV-ANOVA) (%)",
                            "Value": "—",
                            "CI (95%)": "—" if use_cv_anova_flag else "Not computed (enable “Estimate CVI with CV-ANOVA”).",
                        })

                    # remaining rows
                    rows += [
                        {
                            "Parameter": "Between-subject CV (%)",
                            "Value": f"{r.cv_G:.2f} %",
                            "CI (95%)": _fmt_ci(r.ci_cv_G, " %"),
                        },
                        {
                            "Parameter": "Reference change value (95%) (lognormal)",
                            "Value": f"−{r.rcv_95_down:.2f} % / +{r.rcv_95_up:.2f} %",
                            "CI (95%)": "",
                        },
                    ]

                    return pd.DataFrame(rows)



                def _render_key_metrics_table(label: str, r: BVResult, unit: str, use_cv_anova_flag: bool):
                    st.markdown(f"### {label}")

                    # Build the two tables
                    design_df = _build_study_design_table(r)
                    core_df   = _build_core_estimates_table(r, unit, use_cv_anova_flag)

                    # Common style (matches your current header theme)
                    hdr_style = [{
                        "selector": "th",
                        "props": [
                            ("background-color", "#d0e7f5"),
                            ("color", "#1e3a5f"),
                            ("font-size", "0.9rem"),
                        ],
                    }]

                    st.markdown("#### Study design")
                    st.dataframe(
                        design_df.style
                            .set_table_styles(hdr_style)
                            .set_properties(**{"font-size": "0.9rem"})
                            .set_properties(subset=["Parameter"], **{"font-weight": "700"}),
                        hide_index=True,
                        use_container_width=True,
                    )

                    st.markdown("#### Core estimates")
                    st.dataframe(
                        core_df.style
                            .set_table_styles(hdr_style)
                            .set_properties(**{"font-size": "0.9rem"})
                            .set_properties(subset=["Parameter"], **{"font-weight": "700"}),
                        hide_index=True,
                        use_container_width=True,
                    )

                    # Optional note for unbalanced runs
                    try:
                        if (not float(r.S).is_integer()) or (not float(r.R).is_integer()):
                            st.caption("Note: Unbalanced design → S̄ and/or R̄ are reported as mean values.")
                    except Exception:
                        pass



                def _render_key_metrics_section(r: BVResult, unit: str, use_cv_anova_flag: bool):
                    # Build the metrics list (same as your overall)
                    metrics = [
                        ("Mean (CI)",
                        f"{r.grand_mean:.2f}" + (f" {unit}" if unit else ""),
                        f"({r.ci_mean[0]:.2f}–{r.ci_mean[1]:.2f})" + (f" {unit}" if unit else "")),
                        ("CV<sub>A</sub> (CI%)",
                        f"{r.cv_A:.2f}",
                        f"({r.ci_cv_A[0]:.2f}–{r.ci_cv_A[1]:.2f})"),
                        ("CV<sub>I</sub> <span style='font-size:0.7em;'>(based on standard ANOVA)</span> (CI%)",
                        f"{r.cv_I:.2f}",
                        f"({r.ci_cv_I[0]:.2f}–{r.ci_cv_I[1]:.2f})"),
                        ("CV<sub>G</sub> (CI%)",
                        f"{r.cv_G:.2f}",
                        f"({r.ci_cv_G[0]:.2f}–{r.ci_cv_G[1]:.2f})"),
                        ("RCV (95%) (lognormal)",
                        f"−{r.rcv_95_down:.2f}% / +{r.rcv_95_up:.2f}%",
                        None),
                    ]

                    if r.cv_I_cv_anova is not None and r.ci_cv_I_cv_anova is not None:
                        metrics.append((
                            "CV<sub>I</sub> <span style='font-size:0.75em;'>(based on CV-ANOVA)</span> (CI%)",
                            f"{r.cv_I_cv_anova:.2f}",
                            f"({r.ci_cv_I_cv_anova[0]:.2f}–{r.ci_cv_I_cv_anova[1]:.2f})"
                        ))

                    first_row  = metrics[:-3]
                    second_row = metrics[-3:]

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

                    _render_row(first_row)
                    if second_row:
                        st.write(" "); st.write(" ")
                        _render_row(second_row)

                def _render_summary_section(r: BVResult, use_cv_anova_flag: bool | None = None):
                    if use_cv_anova_flag is None:
                        use_cv_anova_flag = st.session_state.get("use_cv_anova", False)
                    # Pick which CVI to display in summary
                    if use_cv_anova_flag and (r.cv_I_cv_anova is not None) and (r.ci_cv_I_cv_anova is not None):
                        cvi_val = float(r.cv_I_cv_anova)
                        cvi_ci  = r.ci_cv_I_cv_anova
                        cvi_lbl = "Within-subject (CVI, CV-ANOVA)"
                    else:
                        cvi_val = float(r.cv_I)
                        cvi_ci  = r.ci_cv_I
                        cvi_lbl = "Within-subject (CVI)"

                    var_tbl = pd.DataFrame({
                        "Variations": ["Analytical", cvi_lbl, "Between-subject"],
                        "Variance":  [r.var_A,    r.var_WP,      r.var_BP],
                        "CV %":      [r.cv_A,     cvi_val,       r.cv_G],
                        "CI (95%)":  [f"{r.ci_cv_A[0]:.2f}–{r.ci_cv_A[1]:.2f}",
                                    f"{cvi_ci[0]:.2f}–{cvi_ci[1]:.2f}",
                                    f"{r.ci_cv_G[0]:.2f}–{r.ci_cv_G[1]:.2f}"],
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



                def _render_aps_section(r: BVResult, use_cv_anova_flag: bool):
                    if use_cv_anova_flag and (r.cv_I_cv_anova is not None) and (r.ci_cv_I_cv_anova is not None):
                        cvi_for_aps = float(r.cv_I_cv_anova)
                        cvi_source  = "CV-ANOVA"
                        cvi_ci      = r.ci_cv_I_cv_anova
                    else:
                        cvi_for_aps = float(r.cv_I)
                        cvi_source  = "standard ANOVA"
                        cvi_ci      = r.ci_cv_I

                    aps_df = build_aps_table(cvi_pct=cvi_for_aps, cvg_pct=r.cv_G)

                    # --- APS TABLE (HTML header so CV subscript A renders) ---
                    aps_df_disp = aps_df.rename(columns={"CVa": "CV<sub>A</sub>"})

                    aps_html = (
                        aps_df_disp.style
                            .format({"CV<sub>A</sub>": "{:.1f}", "Bias": "{:.1f}", "MAU (k=2)": "{:.1f}"})
                            .set_table_styles([
                                {
                                    "selector": "th",
                                    "props": [
                                        ("background-color", "#d0e7f5"),
                                        ("color", "#1e3a5f"),
                                        ("font-weight", "700"),
                                    ],
                                },
                                {"selector": "td", "props": [("text-align", "right")]},
                                {"selector": "td:first-child", "props": [("text-align", "left")]},
                            ])
                            .hide(axis="index")                      # ✅ hide index column
                            .to_html(index=False, escape=False)      # ✅ do not render index
                    )

                    st.markdown(aps_html, unsafe_allow_html=True)


                    st.caption(
                        f"CVI used for APS: **{cvi_source}** → {cvi_for_aps:.2f}% "
                        f"(95% CI {cvi_ci[0]:.2f}–{cvi_ci[1]:.2f}%); "
                        f"MAU: Maximum Allowable Expanded Uncertainty."
                    )

                def _render_by_group(groups: dict[str, tuple[BVResult, pd.DataFrame]], render_one):
                    """
                    If only Overall exists → render once.
                    If Male/Female exist → render in tabs.
                    render_one(label, result, raw_df)
                    """
                    if len(groups) == 1:
                        (label, (r, raw_df)), = groups.items()
                        render_one(label, r, raw_df)
                        return

                    tabs = st.tabs(list(groups.keys()))
                    for tab, (label, (r, raw_df)) in zip(tabs, groups.items()):
                        with tab:
                            render_one(label, r, raw_df)


                # ─────────────────────────────────────────────────────────────
                # Reusable render: Key metrics / Summary / APS / BIVAC
                # ─────────────────────────────────────────────────────────────



                # Key metrics — big bold labels, smaller numbers
                #  Key metrics – two-row layout
                # ────────────────────────────────────────────────────────────
                # st.subheader("Key metrics (cards)")
                # _render_by_group(groups, lambda label, r, raw_df: _render_key_metrics_section(r, unit, use_cv_anova_flag))

                st.write(" ")
                st.subheader("Key metrics (table)")
                _render_by_group(groups, lambda label, r, raw_df: _render_key_metrics_table(label, r, unit, use_cv_anova_flag))



                # Summary table of CVs and 95% CIs
                # st.write(" ")
                # st.write(" ")
                # st.subheader("Summary of variation metrics")

                # _render_by_group(groups, lambda label, r, raw_df: _render_summary_section(r, use_cv_anova_flag))


                # -------------------- APS TABLE --------------------
                st.divider()
                st.write(" ")
                st.subheader("Analytical Performance Specifications (APS)")

                _render_by_group(groups, lambda label, r, raw_df: _render_aps_section(r, use_cv_anova_flag))

                # ---------------------------------------------------



                # — Per-subject mean ± range plot ————————————————————————————————
                # Build only the final dataset unconditionally (balanced/unbalanced per sidebar)
                final_df = res.clean_df
                if final_df is None or final_df.empty:
                    st.error("Cleaned dataset is empty after QC.")
                    st.stop()
                n_raw  = len(df_for_calc)
                n_kept = len(final_df)
                
                st.subheader("Per-subject distribution (overall)")
                st.plotly_chart(plot_subject_ranges(final_df), use_container_width=True)

                if st.session_state.get("gender_based", False) and "df_gender" in locals() and "Gender" in df_gender.columns:
                     st.subheader("Gender-stratified distribution")
                     # Re-map gender to the final cleaned dataframe using the calculation-time gender df
                     gender_map = df_gender.set_index("Subject")["Gender"].to_dict()
                     st.plotly_chart(plot_gender_stratified_ranges(final_df, gender_map), use_container_width=True)
                elif st.session_state.get("gender_based", False):
                     st.warning("Gender plot could not be generated: 'Gender' column missing in processed data.")

                # ⬇️ NEW: run/show population trend ONLY if the sidebar switch is ON
                if st.session_state["preproc_flags"].get("pop_drift", False):
                    # Build the pre-balance dataset only now (so we don't pay the cost unless needed)
                    prebal_df, _ = _preprocess_bv_dataframe(
                        df_for_calc,
                        flags=st.session_state["preproc_flags"],
                        enforce_balance=False
                    )

                    st.subheader("Population time trend (pre-balance)")
                    try:
                        # (Optional micro-optimization: cache this call or reuse a value stored in session_state)
                        pop = estimate_population_drift(prebal_df, alpha=0.05, collapse_replicates=True)
                        unit = st.session_state.get("result_unit", "").strip()
                        unit_txt = f" {unit}" if unit else ""
                        st.markdown(
                            f"**Slope per sample:** {pop['slope']:+.4g}{unit_txt}/sample "
                            f"(95% CI {pop['ci_low']:.4g} to {pop['ci_high']:.4g}, "
                            f"p = {pop['p']:.3g}, df = {pop['df']})."
                        )
                        st.plotly_chart(plot_population_trend(prebal_df), use_container_width=True)
                    except Exception as e:
                        st.info(f"Population drift not computed: {e}")
                else:
                    st.caption("Enable **Population drift (QI 7)** in the sidebar to compute and display the pooled time trend.")




                # ---------------------------------------------------------------------------
                #  📋  BIVAC CHECKLIST  (QI 6 → QI 14)
                #      – fills Comment / Details and skips QI 1‑5
                # ---------------------------------------------------------------------------
                import re
                from collections import defaultdict

                # ── render the checklist & XLSX download ────────────────────────────────────
                st.subheader("BIVAC Checklist (QI 1 → QI 14)")

                def _render_one_bivac(label: str, r: BVResult, raw_df: pd.DataFrame):
                    final_df_x = r.clean_df
                    if final_df_x is None or final_df_x.empty:
                        st.error(f"{label}: cleaned dataset is empty after QC.")
                        return

                    n_raw_x = len(raw_df)
                    n_kept_x = len(final_df_x)

                    unit = st.session_state.get("result_unit", "").strip()

                    bivac_df = _build_qi_checklist(
                        r.preprocess_log,
                        r,
                        flags=flags,
                        qi_manual=qi_manual,
                        n_raw=n_raw_x,
                        n_kept=n_kept_x,
                        unit=unit,   # ✅ keeps units showing
                    )

                    st.dataframe(
                        bivac_df.style
                            .apply(
                                lambda s: ["background:#e8f9f0" if g in ("A", "B") else "background:#fdecea" for g in s],
                                axis=1, subset=["Grade"]
                            )
                            .apply(_style_details_series, subset=["Details"])
                            .set_properties(**{"font-size": "0.85rem"}),
                        use_container_width=True,
                        hide_index=True,
                    )

                    st.markdown(f"**Overall BIVAC grade: {bivac_df.attrs['overall_grade']}**")



                _render_by_group(groups, lambda label, r, raw_df: _render_one_bivac(label, r, raw_df))


                # — Outlier / normality audit trail —
                with st.expander("Pre-processing log"):
                    if st.session_state.get("gender_based", False) and (res_male is not None) and (res_female is not None):
                        t0, t1, t2 = st.tabs(["Overall", "Male", "Female"])
                        with t0:
                            st.markdown(_render_log_html(res.preprocess_log, alpha=0.05), unsafe_allow_html=True)
                        with t1:
                            st.markdown(_render_log_html(res_male.preprocess_log, alpha=0.05), unsafe_allow_html=True)
                        with t2:
                            st.markdown(_render_log_html(res_female.preprocess_log, alpha=0.05), unsafe_allow_html=True)
                    else:
                        st.markdown(_render_log_html(res.preprocess_log, alpha=0.05), unsafe_allow_html=True)

                # ✅ Download ALL checklists in one workbook (multi-sheet)
                with st.expander("⇢ Download BIVAC checklist", expanded=False):
                    import importlib
                    engine = "xlsxwriter" if importlib.util.find_spec("xlsxwriter") else "openpyxl"

                    all_xio = io.BytesIO()
                    with pd.ExcelWriter(all_xio, engine=engine) as xl:
                        for label, (r, raw_df) in groups.items():
                            try:
                                biv = build_bivac_df_for_group(label, r, raw_df)
                                # Excel sheet names max 31 chars
                                sheet = f"BIVAC_{label}"[:31]
                                biv.to_excel(xl, index=False, sheet_name=sheet)
                            except Exception as e:
                                # still write a small sheet explaining why it's missing
                                sheet = f"BIVAC_{label}"[:31]
                                pd.DataFrame({"Error": [str(e)]}).to_excel(xl, index=False, sheet_name=sheet)

                    st.download_button(
                        label="⬇️ Download BIVAC checklist (XLSX)",
                        data=all_xio.getvalue(),
                        file_name="BIVAC_QI_1-14_AllGroups.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        key="dl_bivac_allgroups",
                    )


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

