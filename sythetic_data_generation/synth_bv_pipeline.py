import numpy as np
import pandas as pd


# -------- Targets from your report --------
MU = 0.91              # mean
CV_A = 2.1/ 100      # analytical
CV_I = 4.4 / 100      # within-subject (standard ANOVA)
CV_G = 16.2 / 100     # between-subject

SD_A = MU * CV_A
SD_I = MU * CV_I
SD_G = MU * CV_G

# -------- Dataset shape (960 rows: 20×24×2) --------
N_SUBJECTS = 48
N_SAMPLES_PER_SUBJECT = 10
N_REPS = 2

rng = np.random.default_rng(7)  # reproducible

rows = []
for i in range(1, N_SUBJECTS + 1):
    b_i = rng.normal(0, SD_G)  # between-subject
    for s in range(1, N_SAMPLES_PER_SUBJECT + 1):
        w_is = rng.normal(0, SD_I)  # within-subject
        for r in range(1, N_REPS + 1):
            e_isr = rng.normal(0, SD_A)  # analytical replicate noise
            y = MU + b_i + w_is + e_isr
            rows.append((i, s, r, y))

df = pd.DataFrame(rows, columns=["Subject", "Sample", "Replicate", "Result"])

# Center 
df["Result"] = df["Result"] - (df["Result"].mean() - MU)

# Rounding
df["Result"] = df["Result"].round(2)

# Save
out_path = "crea_simulated_results.xlsx"
with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
    df.to_excel(writer, sheet_name="Sayfa1", index=False)

print("Saved:", out_path)
