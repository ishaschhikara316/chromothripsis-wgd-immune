"""
A1 (dissociation) + A2 (two-channel causal mediation).

A1: mutually adjust chromothripsis & WGD for immune-cold (cytolytic).
A2: decompose WGD -> cytolytic into
      sensor channel  M1 = cGAS-STING sensor activity (ISG + STING1 + CGAS)  [higher = ON]
      APM channel     M2 = antigen presentation (HLA/B2M/TAP/NLRC5/PSMB)     [higher = intact]
    Regression-based interventional decomposition (VanderWeele): indirect via each
    mediator = a-path * b-path; direct = residual. Bootstrap 95% CIs; E-values.
    Confounders: tumour purity, proliferation; tumour type via within-type z-scoring.
Outcome standardised within type => coefficients are ~SD units (Cohen's d).
"""
import warnings, numpy as np, pandas as pd
import statsmodels.formula.api as smf
from pathlib import Path
warnings.filterwarnings("ignore")
rng = np.random.default_rng(20260712)
ROOT = Path(__file__).resolve().parents[1]

chromo = pd.read_csv(ROOT/"data"/"chromo_per_donor.csv").dropna(subset=["barcode"])
def load(fp):
    e = pd.read_csv(ROOT/"data"/fp, sep="\t", index_col=0).T
    e = e[e.index.str.endswith("-01")]; e["barcode"] = e.index.str[:12]
    return e[~e["barcode"].duplicated()]
expr = load("panel_expr.tsv").merge(load("sensor_expr.tsv"), on="barcode", suffixes=("","_s")) \
                             .merge(load("immune_hypoxia_expr.tsv"), on="barcode", suffixes=("","_i"))
df = chromo.merge(expr, on="barcode", how="inner")

ISG=["ISG15","MX1","OAS1","IFIT1","IFIT3","IRF7","DDX58","RSAD2","IFI44","IFI6","USP18","HERC5","STAT1"]
SENSOR = ISG + ["TMEM173","C6orf150"]                       # ISG + STING1 + CGAS
APM=["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","NLRC5","PSMB8","PSMB9"]
PROLIF=["MKI67","PCNA","TOP2A","CCNB1"]
def zmean(frame, genes):
    genes=[g for g in genes if g in frame.columns]; o=pd.DataFrame(index=frame.index)
    for g in genes: o[g]=frame.groupby("histo")[g].transform(lambda v:(v-v.mean())/v.std(ddof=0) if v.std(ddof=0) else 0.0)
    return o.mean(axis=1)
df["CYT"]=zmean(df,["GZMA","PRF1"]); df["M_sensor"]=zmean(df,SENSOR)
df["M_apm"]=zmean(df,APM); df["prolif"]=zmean(df,PROLIF)
df["WGD"]=(df.ploidy>2.5).astype(int); df["chromo"]=df.chromo_high.astype(int)
d=df.dropna(subset=["CYT","M_sensor","M_apm","WGD","chromo","purity","prolif"]).copy().reset_index(drop=True)
print(f"n={len(d)}  WGD+={int(d.WGD.sum())}  chromo+={int(d.chromo.sum())}  types={d.histo.nunique()}\n")

def evalue(est):                     # VanderWeele-Ding E-value for a standardised effect
    rr=np.exp(0.91*abs(est));        # approx RR from Cohen's d
    return rr+np.sqrt(rr*(rr-1))

# ---------- A1: dissociation ----------
print("="*68,"\nA1 — DISSOCIATION: who drives immune-cold, chromothripsis or WGD?\n"+"="*68)
m=smf.ols("CYT ~ chromo + WGD + purity + prolif", data=d).fit()
for v in ["chromo","WGD"]:
    print(f"  {v:7s}: b={m.params[v]:+.3f}  p={m.pvalues[v]:.2e}  E-value={evalue(m.params[v]):.2f}")
print("  -> interpretation: larger |b| / smaller p = dominant driver")

# ---------- A2: two-channel mediation for WGD ----------
def mediation(data, exposure):
    C="purity + prolif"
    tot=smf.ols(f"CYT ~ {exposure} + {C}", data=data).fit().params[exposure]
    a_s=smf.ols(f"M_sensor ~ {exposure} + {C}", data=data).fit().params[exposure]
    a_a=smf.ols(f"M_apm ~ {exposure} + {C}", data=data).fit().params[exposure]
    out=smf.ols(f"CYT ~ {exposure} + M_sensor + M_apm + {C}", data=data).fit()
    b_s, b_a = out.params["M_sensor"], out.params["M_apm"]
    direct=out.params[exposure]
    ind_s, ind_a = a_s*b_s, a_a*b_a
    return dict(total=tot, direct=direct, ind_sensor=ind_s, ind_apm=ind_a,
                pm_sensor=ind_s/tot if tot else np.nan, pm_apm=ind_a/tot if tot else np.nan)

def boot(exposure, B=2000):
    keys=["total","direct","ind_sensor","ind_apm","pm_sensor","pm_apm"]
    acc={k:[] for k in keys}
    for _ in range(B):
        bs=d.sample(len(d), replace=True, random_state=int(rng.integers(1e9)))
        try:
            r=mediation(bs, exposure)
            for k in keys: acc[k].append(r[k])
        except Exception: pass
    return {k:(np.percentile(v,2.5), np.percentile(v,97.5)) for k,v in acc.items()}

for exposure in ["WGD","chromo"]:
    print("\n"+"="*68,f"\nA2 — two-channel mediation, exposure = {exposure}\n"+"="*68)
    pt=mediation(d, exposure); ci=boot(exposure)
    def line(name,key,pm=None):
        lo,hi=ci[key]; extra=""
        if pm is not None:
            plo,phi=ci[pm]; extra=f"   proportion mediated={pt[pm]*100:+.0f}% [{plo*100:+.0f},{phi*100:+.0f}]"
        print(f"  {name:22s} {pt[key]:+.3f} [{lo:+.3f},{hi:+.3f}]{extra}")
    line("total effect","total")
    line("direct effect","direct")
    line("indirect via SENSOR","ind_sensor","pm_sensor")
    line("indirect via APM","ind_apm","pm_apm")
    print(f"  E-value(total)={evalue(pt['total']):.2f}  E-value(sensor)={evalue(pt['ind_sensor']):.2f}  E-value(APM)={evalue(pt['ind_apm']):.2f}")

d.to_csv(ROOT/"results"/"analysis_dataset.csv", index=False)
print("\nwrote results/analysis_dataset.csv")
