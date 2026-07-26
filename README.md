# chromothripsis-wgd-immune

**Within TCGA, whole-genome doubling — not the chromothripsis it travels with — tracks the immune-cold phenotype of complex-genome tumours. But the association DOES NOT REPLICATE in independent non-TCGA cohorts.**

> ⚠️ **Replication failed (A5).** The WGD effect is −0.339 (p=6.6e-7) in TCGA but −0.069 (p=0.56) in non-TCGA PCAWG solid tumours, and this is *not* explained by tumour-type mix, lymphoid contamination, data quality or power (all tested — see `results/RESULTS_CORE.md`). The finding is cohort-dependent and should not be presented as established biology.

A computational, public-data, CPU-only causal-mediation project. Built on **Isidro Cortés-Ciriano's** own published chromothripsis calls (PCAWG ShatterSeek, Cortés-Ciriano et al. *Nat Genet* 2020, PMID 32025003) — his ERC programme *BrokenChromosomes* studies complex rearrangements → immune evasion.

## Headline findings (n = 689 TCGA tumours, 23 types)

1. **Dissociation (the novel core).** Chromothripsis and WGD have never been mutually adjusted. When they are, WGD dominates (b = −0.339, p = 6.6 × 10⁻⁷, E-value 2.06) while chromothripsis retains only a marginal effect (b = −0.144, p = 0.038) that **falls below significance once TMB is adjusted** (p = 0.053). The reported chromothripsis–immune association is largely WGD acting through it.

2. **Not a neoantigen-deficit effect.** WGD tumours carry *more* mutations (b = +0.235, p = 1.1 × 10⁻⁴) yet are colder — so tumour mutational burden cannot explain the phenotype, and adjusting for it *strengthens* the WGD effect.

3. **Antigen presentation, not the innate sensor.** Two-channel decomposition of WGD → immune-cold: antigen-presentation loss mediates **28% [9, 47]**; the **cGAS-STING sensor channel is null (1% [−3, 6])**. The project's original "sensor suppression is the novel channel" hypothesis was **refuted** by formal mediation, and is reported as such. The APM result confirms an independently reported 2026 finding (Foidart et al., *Cancer Cell*, PRC2/MHC-I); the contribution here is causal quantification, not a new mechanism.

4. **~71% of the WGD effect is direct/unexplained.** A mediator hunt found no clean novel channel — the large apparent mediators (Treg, myeloid) are co-linear infiltration markers, and the one directional axis (stroma/TGF-β) does not mediate WGD. The residual reads as global immune desertification that bulk RNA cannot localise.

See `results/RESULTS_CORE.md` for full tables and honest limitations.

## Aims

- **A1 — Dissociation.** Mutually adjust chromothripsis and WGD for immune-cold. ✅ done
- **A2 — Two-channel mediation.** cGAS-STING sensor vs antigen presentation. ✅ done (sensor refuted, APM confirmed)
- **A2b — Residual hunt.** Chase the unexplained direct effect. ✅ done (no clean novel mediator)
- **A3 — TMB deconfounding.** ✅ done (dissociation survives and sharpens)
- **A4 — Figures.** ✅ done
- **A5 — Replication.** PCAWG-direct RNA / ICGC. ✅ done — **does not replicate**

## Data (public, CPU-only)

Cortés-Ciriano 2020 chromothripsis calls (*Nat Genet* supp. Table S1) → TCGA subset via `tcga_donor_uuid`; TCGA RNA (Xena PanCanAtlas); WGD/purity/ploidy/PGA from the chromothripsis supplement (WGD = ploidy > 2.5); TMB from GDC PanImmune mutation-load. All scores z-standardised **within tumour type**, so coefficients are ~SD units.

## Reproduce

```bash
python src/a1_a2_mediation.py     # A1 dissociation + A2 two-channel mediation
python src/a2b_mediator_hunt.py   # A2b residual mediator hunt
python src/a3_tmb_deconfound.py   # A3 TMB deconfounding
python src/a4_figures.py          # A4 figures
```

## What this sets up

The unresolved ~71% direct effect is precisely a **single-cell / spatial** question — the methodology of the group whose chromothripsis calls this is built on. That is the concrete collaboration hook, rather than a generic ask.

*Status: the headline does not replicate. Not suitable for submission as a positive finding; see RESULTS_CORE.md for the full replication analysis.*
