"""
A3 — TMB deconfounding (the outstanding TODO from the A1/A2 build).

Tumour mutational burden is the obvious alternative explanation for an immune-cold
phenotype: fewer mutations -> fewer neoantigens -> colder tumour. WGD and chromothripsis
both co-vary with mutational processes, so the A1 dissociation and the A2 decomposition
are only trustworthy if they survive adjustment for TMB. The original build could not
test this (the TMB download failed); the file is now available locally.

This script:
  1. Reproduces A1 (mutual adjustment of chromothripsis vs WGD) from the saved
     analysis_dataset.csv, to confirm the published numbers still hold.
  2. Joins per-patient TMB (GDC PanImmune mutation-load, "Non-silent per Mb").
  3. Re-runs A1 and the A2 two-channel mediation with log-TMB added as a confounder.

Outputs results/a3_tmb_deconfound.tsv and prints a before/after comparison.
"""
import warnings, numpy as np, pandas as pd, statsmodels.formula.api as smf
from pathlib import Path
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
TMB_SRC = ROOT.parent / "ecdna-immune-mediation" / "data" / "raw" / "mutation-load_updated.txt"
rng = np.random.default_rng(20260719)

d = pd.read_csv(ROOT / "results" / "analysis_dataset.csv")
print(f"loaded n={len(d)}  WGD+={int(d.WGD.sum())}  chromo+={int(d.chromo.sum())}  types={d.histo.nunique()}")


def evalue(est):
    rr = np.exp(0.91 * abs(est))
    return rr + np.sqrt(rr * (rr - 1))


def a1(data, confs):
    m = smf.ols(f"CYT ~ chromo + WGD + {confs}", data=data).fit()
    return {v: (m.params[v], m.pvalues[v], evalue(m.params[v])) for v in ("chromo", "WGD")}


def mediation(data, exposure, confs):
    tot = smf.ols(f"CYT ~ {exposure} + {confs}", data=data).fit().params[exposure]
    a_s = smf.ols(f"M_sensor ~ {exposure} + {confs}", data=data).fit().params[exposure]
    a_a = smf.ols(f"M_apm ~ {exposure} + {confs}", data=data).fit().params[exposure]
    out = smf.ols(f"CYT ~ {exposure} + M_sensor + M_apm + {confs}", data=data).fit()
    ind_s, ind_a = a_s * out.params["M_sensor"], a_a * out.params["M_apm"]
    return dict(total=tot, direct=out.params[exposure], ind_sensor=ind_s, ind_apm=ind_a,
                pm_sensor=ind_s / tot if tot else np.nan, pm_apm=ind_a / tot if tot else np.nan)


def boot(data, exposure, confs, B=2000):
    keys = ["total", "direct", "ind_sensor", "ind_apm", "pm_sensor", "pm_apm"]
    acc = {k: [] for k in keys}
    idx = np.arange(len(data))
    for _ in range(B):
        bs = data.iloc[rng.choice(idx, len(idx), replace=True)]
        try:
            r = mediation(bs, exposure, confs)
            for k in keys:
                acc[k].append(r[k])
        except Exception:
            pass
    return {k: (np.percentile(v, 2.5), np.percentile(v, 97.5)) for k, v in acc.items()}


# ---------- 1. reproduce A1 (no TMB) ----------
BASE = "purity + prolif"
print("\n" + "=" * 70 + "\nA1 REPRODUCTION (confounders: purity + proliferation)\n" + "=" * 70)
r_base = a1(d, BASE)
for v, (b, p, e) in r_base.items():
    print(f"  {v:7s} b={b:+.3f}  p={p:.2e}  E-value={e:.2f}")

# ---------- 2. join TMB ----------
tmb = pd.read_csv(TMB_SRC, sep="\t")
tmb = tmb.rename(columns={"Patient_ID": "barcode", "Non-silent per Mb": "tmb"})[["barcode", "tmb"]]
tmb["tmb"] = pd.to_numeric(tmb["tmb"], errors="coerce")
tmb = tmb.dropna().drop_duplicates("barcode")
d2 = d.merge(tmb, on="barcode", how="left")
d2["log_tmb"] = np.log1p(d2["tmb"])
matched = d2["log_tmb"].notna().sum()
print(f"\nTMB joined: {matched}/{len(d2)} tumours ({matched/len(d2):.0%})")
dt = d2.dropna(subset=["log_tmb"]).copy()

# is TMB actually related to the exposures / outcome? (does it deserve adjustment)
print("  sanity — TMB associations:")
for term, formula in [("WGD", "log_tmb ~ WGD + purity"), ("chromo", "log_tmb ~ chromo + purity"),
                      ("CYT", "CYT ~ log_tmb + purity + prolif")]:
    f = smf.ols(formula, data=dt).fit()
    key = term if term != "CYT" else "log_tmb"
    print(f"    {formula:36s} b({key})={f.params[key]:+.3f} p={f.pvalues[key]:.2e}")

# ---------- 3. A1 + A2 with TMB ----------
FULL = "purity + prolif + log_tmb"
print("\n" + "=" * 70 + "\nA1 WITH TMB ADJUSTMENT\n" + "=" * 70)
r_tmb = a1(dt, FULL)
rows = []
for v in ("chromo", "WGD"):
    b0, p0, e0 = a1(dt, BASE)[v]          # same n, no TMB -> fair comparison
    b1, p1, e1 = r_tmb[v]
    print(f"  {v:7s} before b={b0:+.3f} p={p0:.2e}   ->   after b={b1:+.3f} p={p1:.2e}  E={e1:.2f}")
    rows.append(dict(aim="A1", term=v, n=len(dt), b_noTMB=b0, p_noTMB=p0, b_TMB=b1, p_TMB=p1, evalue_TMB=e1))

print("\n" + "=" * 70 + "\nA2 TWO-CHANNEL MEDIATION (WGD) WITH TMB\n" + "=" * 70)
pt = mediation(dt, "WGD", FULL)
ci = boot(dt, "WGD", FULL)
for name, key, pmk in [("total effect", "total", None), ("direct effect", "direct", None),
                       ("indirect via SENSOR", "ind_sensor", "pm_sensor"),
                       ("indirect via APM", "ind_apm", "pm_apm")]:
    lo, hi = ci[key]
    extra = ""
    if pmk:
        plo, phi = ci[pmk]
        extra = f"   proportion mediated={pt[pmk]*100:+.0f}% [{plo*100:+.0f},{phi*100:+.0f}]"
    print(f"  {name:22s} {pt[key]:+.3f} [{lo:+.3f},{hi:+.3f}]{extra}")
    rows.append(dict(aim="A2_WGD_TMBadj", term=name, n=len(dt), est=pt[key], lo=lo, hi=hi,
                     pm=pt[pmk] if pmk else np.nan))

pd.DataFrame(rows).to_csv(ROOT / "results" / "a3_tmb_deconfound.tsv", sep="\t", index=False)
print("\nwrote results/a3_tmb_deconfound.tsv")
