"""
Hunt the ~71% direct WGD->immune-cold effect: screen candidate T-cell-exclusion
programs as mediators, alongside APM (known) and the sensor (null).
Single-mediator proportion-mediated (bootstrap CI) + a full multi-mediator model.
"""
import warnings, numpy as np, pandas as pd
import statsmodels.formula.api as smf
from pathlib import Path
warnings.filterwarnings("ignore")
rng=np.random.default_rng(7)
ROOT=Path(__file__).resolve().parents[1]

chromo=pd.read_csv(ROOT/"data"/"chromo_per_donor.csv").dropna(subset=["barcode"])
def load(fp):
    e=pd.read_csv(ROOT/"data"/fp,sep="\t",index_col=0).T; e=e[e.index.str.endswith("-01")]
    e["barcode"]=e.index.str[:12]; return e[~e["barcode"].duplicated()]
expr=load("panel_expr.tsv")
for fp in ["sensor_expr.tsv","immune_hypoxia_expr.tsv","mediator_expr.tsv"]:
    expr=expr.merge(load(fp),on="barcode",suffixes=("","_"+fp[:3]))
df=chromo.merge(expr,on="barcode",how="inner")

gsets={ln.split("\t")[0]:ln.rstrip().split("\t")[1].split(",") for ln in open(ROOT/"data"/"mediator_genesets.txt")}
ISG=["ISG15","MX1","OAS1","IFIT1","IFIT3","IRF7","DDX58","RSAD2","IFI44","IFI6","USP18","HERC5","STAT1"]
progs={"APM":["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","NLRC5","PSMB8","PSMB9"],
       "SENSOR":ISG+["TMEM173","C6orf150"], **gsets}
def zmean(frame,genes):
    genes=[g for g in genes if g in frame.columns]; o=pd.DataFrame(index=frame.index)
    for g in genes: o[g]=frame.groupby("histo")[g].transform(lambda v:(v-v.mean())/v.std(ddof=0) if v.std(ddof=0) else 0.0)
    return o.mean(axis=1)
df["CYT"]=zmean(df,["GZMA","PRF1"]); df["prolif"]=zmean(df,["MKI67","PCNA","TOP2A","CCNB1"])
for p,g in progs.items(): df["M_"+p]=zmean(df,g)
df["WGD"]=(df.ploidy>2.5).astype(int)
mcols=["M_"+p for p in progs]
d=df.dropna(subset=["CYT","WGD","purity","prolif"]+mcols).copy().reset_index(drop=True)
print(f"n={len(d)}  WGD+={int(d.WGD.sum())}\n")

C="purity + prolif"
total=smf.ols(f"CYT ~ WGD + {C}",data=d).fit().params["WGD"]
print(f"total WGD->CYT = {total:+.3f}\n")
print("=== single-mediator screen: proportion of WGD->immune-cold via each program ===")
print(f"{'program':9s} {'a(WGD->M)':>10s} {'b(M->CYT)':>10s} {'%mediated':>10s} {'95% CI':>16s}")
def pm_boot(M,B=1500):
    vals=[]
    for _ in range(B):
        bs=d.sample(len(d),replace=True,random_state=int(rng.integers(1e9)))
        try:
            t=smf.ols(f"CYT ~ WGD + {C}",data=bs).fit().params["WGD"]
            a=smf.ols(f"{M} ~ WGD + {C}",data=bs).fit().params["WGD"]
            b=smf.ols(f"CYT ~ WGD + {M} + {C}",data=bs).fit().params[M]
            vals.append((a*b/t)*100 if t else np.nan)
        except Exception: pass
    return np.nanpercentile(vals,2.5),np.nanpercentile(vals,97.5)
rows=[]
for p in progs:
    M="M_"+p
    a=smf.ols(f"{M} ~ WGD + {C}",data=d).fit().params["WGD"]
    b=smf.ols(f"CYT ~ WGD + {M} + {C}",data=d).fit().params[M]
    pm=(a*b/total)*100
    lo,hi=pm_boot(M)
    rows.append((p,a,b,pm,lo,hi))
for p,a,b,pm,lo,hi in sorted(rows,key=lambda r:-r[3]):
    sig="*" if (lo>0 or hi<0) else " "
    print(f"{p:9s} {a:+10.3f} {b:+10.3f} {pm:+9.0f}% {f'[{lo:+.0f},{hi:+.0f}]':>16s} {sig}")

print("\n=== full multi-mediator model (all programs jointly): direct effect left ===")
full=smf.ols(f"CYT ~ WGD + {' + '.join(mcols)} + {C}",data=d).fit()
print(f"  direct WGD effect with ALL mediators = {full.params['WGD']:+.3f} (p={full.pvalues['WGD']:.1e})")
print(f"  => jointly mediated = {(1-full.params['WGD']/total)*100:.0f}% of total")
print("  mediator b-paths (effect on CYT | others):")
for p in progs:
    M="M_"+p; print(f"    {p:9s} b={full.params[M]:+.3f} p={full.pvalues[M]:.1e}")
