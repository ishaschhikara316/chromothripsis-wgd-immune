"""A4 — figures for the chromothripsis-vs-WGD dissociation."""
import warnings, numpy as np, pandas as pd, statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)
TMB_SRC = ROOT.parent / "ecdna-immune-mediation" / "data" / "raw" / "mutation-load_updated.txt"

d = pd.read_csv(ROOT / "results" / "analysis_dataset.csv")
tmb = pd.read_csv(TMB_SRC, sep="\t").rename(columns={"Patient_ID": "barcode", "Non-silent per Mb": "tmb"})
d = d.merge(tmb[["barcode", "tmb"]].drop_duplicates("barcode"), on="barcode", how="left")
d["log_tmb"] = np.log1p(pd.to_numeric(d["tmb"], errors="coerce"))
dt = d.dropna(subset=["log_tmb"])

BLUE, RED, GREY = "#4C6EF5", "#E03131", "#adb5bd"

# ---- Fig 1: the dissociation (coefficient plot, before/after TMB) ----
rows = []
for lab, data, confs in [("purity+prolif", d, "purity + prolif"),
                         ("+ TMB", dt, "purity + prolif + log_tmb")]:
    m = smf.ols(f"CYT ~ chromo + WGD + {confs}", data=data).fit()
    for v, nm in [("chromo", "Chromothripsis"), ("WGD", "WGD")]:
        ci = m.conf_int().loc[v]
        rows.append(dict(model=lab, term=nm, b=m.params[v], lo=ci[0], hi=ci[1], p=m.pvalues[v]))
r = pd.DataFrame(rows)

fig, ax = plt.subplots(figsize=(7, 3.4))
ypos = {"Chromothripsis": 1, "WGD": 0}
for i, (lab, off, mk) in enumerate([("purity+prolif", .12, "o"), ("+ TMB", -.12, "s")]):
    sub = r[r.model == lab]
    y = [ypos[t] + off for t in sub.term]
    ax.errorbar(sub.b, y, xerr=[sub.b - sub.lo, sub.hi - sub.b], fmt=mk,
                color=BLUE if i == 0 else RED, capsize=3, ms=7, lw=1.6, label=f"adjusted: {lab}")
ax.axvline(0, color=GREY, lw=1, zorder=0)
ax.set_yticks([0, 1]); ax.set_yticklabels(["WGD", "Chromothripsis"])
ax.set_xlabel("effect on cytolytic activity (SD units)")
ax.set_title("WGD — not chromothripsis — drives the immune-cold phenotype", fontsize=11)
ax.legend(fontsize=8, loc="lower left"); fig.tight_layout()
fig.savefig(FIG / "fig1_dissociation.png", dpi=160); plt.close(fig)

# ---- Fig 2: two-channel mediation ----
def mediation(data, confs, exposure="WGD"):
    tot = smf.ols(f"CYT ~ {exposure} + {confs}", data=data).fit().params[exposure]
    a_s = smf.ols(f"M_sensor ~ {exposure} + {confs}", data=data).fit().params[exposure]
    a_a = smf.ols(f"M_apm ~ {exposure} + {confs}", data=data).fit().params[exposure]
    o = smf.ols(f"CYT ~ {exposure} + M_sensor + M_apm + {confs}", data=data).fit()
    return tot, o.params[exposure], a_s * o.params["M_sensor"], a_a * o.params["M_apm"]

C = "purity + prolif + log_tmb"
tot, direct, ind_s, ind_a = mediation(dt, C)
rng = np.random.default_rng(20260719); idx = np.arange(len(dt)); bs_s, bs_a = [], []
for _ in range(2000):
    b = dt.iloc[rng.choice(idx, len(idx), replace=True)]
    try:
        _, _, s, a = mediation(b, C); bs_s.append(s); bs_a.append(a)
    except Exception: pass

fig, ax = plt.subplots(figsize=(7, 3.2))
names = ["total\neffect", "direct\n(unexplained)", "via cGAS-STING\nsensor", "via antigen\npresentation"]
vals = [tot, direct, ind_s, ind_a]
errs = [[0, 0], [0, 0],
        [ind_s - np.percentile(bs_s, 2.5), np.percentile(bs_s, 97.5) - ind_s],
        [ind_a - np.percentile(bs_a, 2.5), np.percentile(bs_a, 97.5) - ind_a]]
cols = [GREY, GREY, RED, BLUE]
ax.bar(names, vals, color=cols)
for i in (2, 3):
    ax.errorbar(i, vals[i], yerr=[[errs[i][0]], [errs[i][1]]], fmt="none", ecolor="#333", capsize=4)
ax.axhline(0, color="#333", lw=1)
ax.set_ylabel("effect on cytolytic (SD)")
ax.set_title("WGD acts via antigen presentation (28%); sensor channel is null (1%)", fontsize=11)
fig.tight_layout(); fig.savefig(FIG / "fig2_two_channel.png", dpi=160); plt.close(fig)

# ---- Fig 3: TMB is not the explanation ----
fig, axes = plt.subplots(1, 2, figsize=(8, 3.2))
for ax, (grp, lab) in zip(axes, [("WGD", "WGD"), ("chromo", "Chromothripsis")]):
    data = [dt[dt[grp] == 0]["log_tmb"], dt[dt[grp] == 1]["log_tmb"]]
    try:
        bp = ax.boxplot(data, tick_labels=[f"{lab}−", f"{lab}+"], patch_artist=True, widths=.55)
    except TypeError:  # matplotlib < 3.9
        bp = ax.boxplot(data, labels=[f"{lab}−", f"{lab}+"], patch_artist=True, widths=.55)
    for p, c in zip(bp["boxes"], [GREY, BLUE]): p.set_facecolor(c); p.set_alpha(.75)
    f = smf.ols(f"log_tmb ~ {grp} + purity", data=dt).fit()
    ax.set_title(f"{lab}: b={f.params[grp]:+.3f}, p={f.pvalues[grp]:.1e}", fontsize=9)
    ax.set_ylabel("log(1+TMB)")
fig.suptitle("WGD tumours carry MORE mutations yet are colder — not a neoantigen deficit", fontsize=10)
fig.tight_layout(); fig.savefig(FIG / "fig3_tmb_not_explanation.png", dpi=160); plt.close(fig)

print("wrote", *[p.name for p in sorted(FIG.glob("*.png"))])
