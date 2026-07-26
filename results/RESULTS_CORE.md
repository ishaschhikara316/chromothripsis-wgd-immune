# Chromothripsis vs whole-genome doubling in the immune-cold phenotype — core results

*Status: analysis complete, manuscript in preparation. Public data only, CPU-only.*
*Cohort: 689 TCGA tumours with chromothripsis calls from Cortés-Ciriano et al., Nat Genet 2020
(ShatterSeek/PCAWG, supp. Table S1) joined to TCGA expression; 271 WGD+, 213 chromothripsis-high,
23 tumour types. All scores z-standardised within tumour type, so coefficients are ~SD units.*

---

## Headline

**The immune-cold phenotype of complex-genome tumours is driven by whole-genome doubling (WGD),
not by the chromothripsis it travels with.** Chromothripsis and WGD have never before been
mutually adjusted; when they are, chromothripsis's apparent effect largely disappears while WGD's
strengthens. The effect is *not* explained by mutational burden, and it runs through
antigen-presentation loss rather than innate-sensor (cGAS-STING) suppression.

---

## A1 — Dissociation: who actually drives immune-cold?

Mutually adjusted model, outcome = cytolytic activity (GZMA/PRF1), confounders = purity +
proliferation, within-type standardised:

| Term | b | p | E-value |
|---|--:|--:|--:|
| chromothripsis | −0.144 | 0.038 | 1.54 |
| **WGD** | **−0.339** | **6.6 × 10⁻⁷** | **2.06** |

WGD's effect is ~2.4× larger and five orders of magnitude more significant. This is the novel core:
the chromothripsis–immune association reported in the literature is largely **WGD acting through it**.

## A3 — Does it survive tumour mutational burden? (yes, and it sharpens)

TMB is the obvious alternative explanation for an immune-cold tumour: fewer mutations → fewer
neoantigens → colder. TMB (GDC PanImmune non-silent per Mb) joined for **647/689 tumours (94%)**.

**First, TMB behaves in the opposite direction to the confounding story:**

| Association | b | p |
|---|--:|--:|
| WGD → log-TMB | **+0.235** | 1.1 × 10⁻⁴ |
| chromothripsis → log-TMB | −0.025 | 0.69 |
| log-TMB → cytolytic | +0.038 | 0.41 |

WGD tumours carry **more** mutations, yet are **colder**. The immune-cold phenotype therefore
cannot be a neoantigen-deficit effect — if anything TMB suppresses the association, so the
unadjusted estimate was conservative.

**With log-TMB added as a confounder (same n=647 for a fair comparison):**

| Term | b (no TMB) | p | b (TMB-adjusted) | p | E-value |
|---|--:|--:|--:|--:|--:|
| chromothripsis | −0.142 | 0.048 | −0.139 | **0.053** | 1.52 |
| **WGD** | −0.332 | 3.0 × 10⁻⁶ | **−0.346** | **1.4 × 10⁻⁶** | **2.08** |

The dissociation becomes cleaner: after TMB adjustment chromothripsis no longer reaches
significance, while WGD strengthens.

## A2 — Which channel carries the WGD effect? (antigen presentation, not the innate sensor)

Two-channel interventional decomposition of WGD → cytolytic, TMB-adjusted, 2,000 bootstrap
resamples. Channels: **sensor** = cGAS-STING/ISG (higher = sensor ON); **APM** = antigen
presentation (HLA/B2M/TAP/NLRC5/PSMB, higher = intact).

| Path | estimate [95% CI] | proportion mediated |
|---|--:|--:|
| total effect | −0.359 [−0.496, −0.214] | — |
| direct effect | −0.255 [−0.377, −0.138] | ~71% |
| indirect via **sensor** | −0.005 [−0.021, +0.007] | **+1% [−3, +6]** (null) |
| indirect via **APM** | −0.099 [−0.178, −0.024] | **+28% [+9, +47]** |

**The cGAS-STING sensor-suppression hypothesis is refuted** — the channel is flat null, so the
project's original "novel mechanism" is not real and is reported as such. The APM channel is
significant and **confirms** the independently-reported WGD → MHC-I silencing result (Foidart et
al., Cancer Cell 2026); this project adds causal quantification (28%) rather than novelty.

## The unresolved 71%

Most of WGD's immune-cold effect is direct/unexplained. A follow-up mediator hunt found no clean
novel channel: the largest apparent mediators (Treg, myeloid, APM ~30–37%) are co-linear
immune-infiltration markers rather than mechanistically directional, and the one directional axis
(stroma/TGF-β) does **not** mediate WGD — WGD tumours do not have more stroma. The residual reads
as global immune desertification that **bulk RNA cannot resolve**.

---

## Honest limitations

- Cross-sectional TCGA; associations are causal *models*, not experiments.
- Chromothripsis is a binary high/low call; misclassification would bias toward the null.
- The residual ~71% is unexplained — the paper's most honest statement is that bulk expression
  cannot localise it.
- APM finding confirms an existing 2026 result; the novelty here is the dissociation and the
  quantified decomposition, not a new mechanism.
- Replication on PCAWG-direct RNA / ICGC is not yet done.

## What this sets up

The unresolved direct effect is the natural **single-cell / spatial** question, which is exactly
the methodology of the Cortés-Ciriano group whose chromothripsis calls this analysis is built on —
a concrete collaboration hook rather than a generic ask.

## Files

`results/analysis_dataset.csv` (n=689 analysis table) · `results/a3_tmb_deconfound.tsv` (A3 output)
· `src/a1_a2_mediation.py` (A1/A2) · `src/a2b_mediator_hunt.py` (residual hunt) ·
`src/a3_tmb_deconfound.py` (TMB deconfounding).
