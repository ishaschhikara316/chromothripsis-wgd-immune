"""
A5 — Independent replication of the chromothripsis-vs-WGD dissociation in NON-TCGA PCAWG donors.

The discovery cohort (A1-A3) was the TCGA subset of the ShatterSeek chromothripsis calls.
This replicates in the *non-TCGA* (ICGC) donors carrying the SAME chromothripsis calls but a
completely independent expression pipeline (PCAWG tophat/STAR FPKM-UQ, not TCGA Xena).

Note on scope (verified 2026-07-19): PCAWG only ever generated RNA-seq for 7 non-TCGA projects
(LIRI-JP, RECA-EU, MALY-DE, OV-AU, PACA-AU, CLLE-ES, ESAD-UK), so the replication cohort is a few
hundred donors, not the full 1,666 non-TCGA donors with chromothripsis calls. That is a sequencing
-history limit, not an access limit.

Gene symbol -> Ensembl mapping is taken from a local annotated dataset (no hardcoded IDs).
Outputs results/a5_replication.tsv.
"""
import warnings, gzip, numpy as np, pandas as pd, statsmodels.formula.api as smf
from pathlib import Path
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
SSD = Path.home() / "isha-ssd" / "pcawg-replication" / "data" / "raw"
EXPR = SSD / "tophat_star_fpkm_uq.v2_aliquot_gl.tsv.gz"
SHEET = SSD / "pcawg_sample_sheet.tsv"
ATLAS = Path.home() / "isha-ssd" / "pd-sex-dimorphism" / "data" / "raw" / "pd_atlas.h5ad"

ISG = ["ISG15","MX1","OAS1","IFIT1","IFIT3","IRF7","DDX58","RSAD2","IFI44","IFI6","USP18","HERC5","STAT1"]
SENSOR = ISG + ["TMEM173","STING1","CGAS","MB21D1","C6orf150"]
APM = ["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","NLRC5","PSMB8","PSMB9"]
PROLIF = ["MKI67","PCNA","TOP2A","CCNB1"]
CYT = ["GZMA","PRF1"]
WANT = sorted(set(CYT + SENSOR + APM + PROLIF))

# ---------- 1. symbol -> Ensembl from a local annotated dataset ----------
import anndata as ad
v = ad.read_h5ad(ATLAS, backed="r").var
sym2ens = {}
for ens, sym in zip(v.index, v["feature_name"].astype(str)):
    if sym in WANT:
        sym2ens.setdefault(sym, ens.split(".")[0])
print(f"symbol->Ensembl resolved for {len(sym2ens)}/{len(WANT)} target genes")
ens2sym = {e: s for s, e in sym2ens.items()}

# ---------- 2. stream the PCAWG matrix, keep only target genes ----------
rows, header = [], None
with gzip.open(EXPR, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        gid = line.split("\t", 1)[0].split(".")[0]
        if gid in ens2sym:
            rows.append(line.rstrip("\n").split("\t"))
expr = pd.DataFrame(rows, columns=header).set_index("feature")
expr.index = [ens2sym[i.split(".")[0]] for i in expr.index]
expr = expr.astype(float)
expr = np.log2(expr + 1)                      # FPKM-UQ -> log space
print(f"expression panel: {expr.shape[0]} genes x {expr.shape[1]} aliquots")

# ---------- 3. aliquot -> donor ----------
sheet = pd.read_csv(SHEET, sep="\t", low_memory=False)
rna = sheet[sheet["library_strategy"].eq("RNA-Seq")].copy()
spec = rna["dcc_specimen_type"].astype(str).str.lower()
rna = rna[spec.str.contains("tumour|tumor|primary|metasta|recurr", na=False)]
amap = rna.drop_duplicates("aliquot_id").set_index("aliquot_id")["donor_unique_id"].to_dict()

e = expr.T
e["donor_unique_id"] = [amap.get(a) for a in e.index]
e = e.dropna(subset=["donor_unique_id"]).groupby("donor_unique_id").mean()   # collapse to donor
print(f"donors with expression: {len(e)}")

# ---------- 4. join chromothripsis calls, keep NON-TCGA only ----------
chromo = pd.read_csv(ROOT / "data" / "chromo_per_donor.csv")
non_tcga = chromo[chromo["barcode"].isna()].copy()
print(f"non-TCGA donors with chromothripsis calls: {len(non_tcga)}")
d = non_tcga.merge(e, left_on="donor_unique_id", right_index=True, how="inner")
print(f"REPLICATION COHORT n = {len(d)}  |  histologies: {d['histo'].nunique()}")
print(d["histo"].value_counts().to_dict())

# ---------- 5. scores (within-histology z, matching the discovery pipeline) ----------
def zmean(frame, genes):
    genes = [g for g in genes if g in frame.columns]
    o = pd.DataFrame(index=frame.index)
    for g in genes:
        o[g] = frame.groupby("histo")[g].transform(lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) else 0.0)
    return o.mean(axis=1)

d["CYT"] = zmean(d, CYT); d["M_sensor"] = zmean(d, SENSOR)
d["M_apm"] = zmean(d, APM); d["prolif"] = zmean(d, PROLIF)
d["WGD"] = (d["ploidy"] > 2.5).astype(int)
d["chromo"] = d["chromo_high"].astype(int)
d = d.dropna(subset=["CYT", "M_sensor", "M_apm", "WGD", "chromo", "purity", "prolif"])
print(f"analysable n = {len(d)}  WGD+={int(d.WGD.sum())}  chromo+={int(d.chromo.sum())}")

def evalue(est):
    rr = np.exp(0.91 * abs(est)); return rr + np.sqrt(rr * (rr - 1))

# ---------- 6. A1 replication ----------
print("\n" + "=" * 68 + "\nA5 — REPLICATION of A1 dissociation (non-TCGA PCAWG)\n" + "=" * 68)
m = smf.ols("CYT ~ chromo + WGD + purity + prolif", data=d).fit()
out = []
for vname in ("chromo", "WGD"):
    b, p = m.params[vname], m.pvalues[vname]
    ci = m.conf_int().loc[vname]
    print(f"  {vname:7s} b={b:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]  p={p:.3g}  E={evalue(b):.2f}")
    out.append(dict(cohort="PCAWG_nonTCGA", aim="A1", term=vname, n=len(d),
                    b=b, lo=ci[0], hi=ci[1], p=p, evalue=evalue(b)))

# Lineage sensitivity: the discovery cohort was solid tumours. 38% of this replication cohort is
# LYMPHOID (Lymph-BNHL/Lymph-CLL), where "cytolytic activity" (GZMA/PRF1) does not mean the same
# thing — the malignant cell is itself a lymphocyte. Reported alongside the pre-specified analysis.
print("\n  lineage sensitivity (discovery cohort was solid tumours):")
for lab, sub in [("SOLID only", d[~d.histo.str.startswith("Lymph")]),
                 ("LYMPHOID only", d[d.histo.str.startswith("Lymph")])]:
    if len(sub) > 30 and sub.WGD.nunique() > 1:
        mm = smf.ols("CYT ~ chromo + WGD + purity + prolif", data=sub).fit()
        cim = mm.conf_int()
        print(f"    {lab:14s} n={len(sub):3d} WGD+={int(sub.WGD.sum()):3d} | "
              f"WGD b={mm.params['WGD']:+.3f} [{cim.loc['WGD',0]:+.3f},{cim.loc['WGD',1]:+.3f}] "
              f"p={mm.pvalues['WGD']:.3g} | chromo b={mm.params['chromo']:+.3f} p={mm.pvalues['chromo']:.3g}")
        out.append(dict(cohort=f"PCAWG_nonTCGA_{lab.split()[0]}", aim="A1_lineage", term="WGD",
                        n=len(sub), b=mm.params["WGD"], lo=cim.loc["WGD", 0], hi=cim.loc["WGD", 1],
                        p=mm.pvalues["WGD"], evalue=evalue(mm.params["WGD"])))

# per-histology WGD effect (is it consistent?)
print("\n  per-histology WGD effect:")
for h, g in d.groupby("histo"):
    if len(g) >= 30 and g.WGD.nunique() > 1:
        mm = smf.ols("CYT ~ chromo + WGD + purity + prolif", data=g).fit()
        print(f"    {h:20s} n={len(g):3d}  WGD b={mm.params['WGD']:+.3f} p={mm.pvalues['WGD']:.3g}")
        out.append(dict(cohort=h, aim="A1_per_histology", term="WGD", n=len(g),
                        b=mm.params["WGD"], lo=np.nan, hi=np.nan, p=mm.pvalues["WGD"], evalue=np.nan))

# ---------- 7. A2 two-channel mediation in the replication cohort ----------
print("\n" + "=" * 68 + "\nA5 — two-channel mediation (WGD), replication cohort\n" + "=" * 68)
C = "purity + prolif"
tot = smf.ols(f"CYT ~ WGD + {C}", data=d).fit().params["WGD"]
a_s = smf.ols(f"M_sensor ~ WGD + {C}", data=d).fit().params["WGD"]
a_a = smf.ols(f"M_apm ~ WGD + {C}", data=d).fit().params["WGD"]
o = smf.ols(f"CYT ~ WGD + M_sensor + M_apm + {C}", data=d).fit()
ind_s, ind_a = a_s * o.params["M_sensor"], a_a * o.params["M_apm"]
for nm, val in [("total", tot), ("direct", o.params["WGD"]),
                ("indirect SENSOR", ind_s), ("indirect APM", ind_a)]:
    pm = f"  ({val/tot*100:+.0f}% of total)" if tot and nm.startswith("indirect") else ""
    print(f"  {nm:18s} {val:+.3f}{pm}")
    out.append(dict(cohort="PCAWG_nonTCGA", aim="A2", term=nm, n=len(d), b=val,
                    lo=np.nan, hi=np.nan, p=np.nan, evalue=np.nan))

pd.DataFrame(out).to_csv(ROOT / "results" / "a5_replication.tsv", sep="\t", index=False)
print("\nwrote results/a5_replication.tsv")
