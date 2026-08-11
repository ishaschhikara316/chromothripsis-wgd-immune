# Thorough analysis — chromothripsis-wgd-immune and its full lineage

*Critical synthesis across the whole investigation. Complements the results tables in
`results/RESULTS_CORE.md` — this document is the honest "what does it all amount to" assessment.*

---

## Verdict in one paragraph

This was a rigorously executed investigation that reached an honest negative. Within TCGA it
produced a clean, TMB-robust dissociation — the immune-cold phenotype of complex-genome tumours
tracks **whole-genome doubling (WGD), not the chromothripsis it co-occurs with** — but that headline
**does not replicate in an independent non-TCGA cohort**, and every benign explanation for the
non-replication was tested and ruled out. The mechanistic novelty that originally motivated the
project (cGAS-STING **sensor** suppression as a new immune-evasion channel) was **refuted** by formal
mediation; the channel that does carry the effect (antigen presentation) merely **confirms** an
already-published 2026 result. The correct scientific status is therefore: **not a positive flagship
finding.** The genuine, defensible yield is (1) a *cautionary* result about the cohort-dependence of
aneuploidy/WGD–immune associations, (2) a demonstration that the popular "chromothripsis → immune
evasion" association is largely WGD acting through it, and (3) a reusable data resource from the
project's origin (a validated pan-cancer calcification label). None of this required a wet lab; all
of it is public-data, CPU-only.

---

## 1. The full arc — how the project got here

The chromothripsis project was not the starting point. It was the surviving branch of a wider,
disciplined search whose defining feature was **empirical de-risking before commitment**. The arc:

1. **Origin — "calcification around tumour cells" (→ null).** Isha's abstract idea, reframed
   computationally (does an ectopic-mineralization program mediate immune exclusion / outcome
   pan-cancer). A validated calcification label was mined from TCGA pathology reports (TCGA-Reports;
   951 positive / 10%, 140 psammoma), and it tracks known biology. But the core hypotheses were
   **null / heterogeneous**: the osteomimicry signature is a weak proxy for measured calcification
   (CV-AUROC ~0.58; two of Yang 2020's three genes fail); the immune association is null and
   *flips sign* by cancer type (ovarian cold p=0.016, thyroid hot p=0.012); survival is null after
   age adjustment and also flips (uterine worse HR 2.40). Calcification is cancer-type-specific,
   not a unified pan-cancer phenomenon.

2. **Phenotype screen (→ strong survival, null immune).** Reusing the report-mining machinery,
   necrosis / LVI / PNI / sarcomatoid / signet-ring were screened. All carry **strong, clean
   pan-cancer survival signals** (HR 1.25–2.01, p ≤ 1e-7) but the immune association is null across
   the board (~0.50), including the TILs positive control at only 0.53 — i.e. bulk-RNA immune signal
   for report-mined phenotypes is weak, and the survival signals are textbook (aggressive histology →
   death), not novel.

3. **"Signatures = free histology?" (→ refuted).** A proxyaudit-style test of whether molecular
   prognostic signatures are redundant with pathology-report histology. The a-priori strongest case
   (necrosis vs a hypoxia signature) showed necrosis explains only **~5%** of the signature's
   prognostic value (proliferation ~3%). Signatures add independent value; the provocative premise
   was wrong.

4. **PI-anchored idea generation (→ the chromothripsis project).** Instead of guessing in a vacuum,
   idea generation was anchored to the target PIs' work. The primary anchor (Cortés-Ciriano, ERC
   *BrokenChromosomes* = complex rearrangements → immune evasion) pointed directly at chromothripsis,
   and his **own public data** (PCAWG ShatterSeek calls, *Nat Genet* 2020) supplied the exposure.
   This was the first idea in the whole search to survive an empirical de-risk.

**Lesson embedded in the arc:** the cheap up-front de-risking (a few CPU-hours each) correctly killed
three ideas before weeks were sunk into them. That discipline is the transferable asset — more than
any single result.

---

## 2. The chromothripsis-WGD results (A1–A5), critically read

Full tables live in `results/RESULTS_CORE.md`. The critical reading:

- **A1 — dissociation (the novel core, TCGA).** Mutually adjusted (never done before in the
  literature): WGD b=−0.339, p=6.6e-7, E=2.06; chromothripsis b=−0.144, p=0.038. WGD is ~2.4× larger
  and five orders of magnitude more significant. The reported chromothripsis–immune link is largely
  **WGD acting through it**. *This is the genuinely novel analytical contribution.*

- **A3 — TMB (sharpens the dissociation, and a good mechanistic line).** With log-TMB added
  (n=647/689): WGD strengthens to −0.346 (p=1.4e-6, E=2.08); chromothripsis drops to p=0.053 —
  **no independent effect once TMB is controlled.** The reusable line: WGD tumours carry *more*
  mutations (b=+0.235, p=1.1e-4) yet are *colder*, and TMB→cytolytic is null — so the phenotype is
  **not a neoantigen-deficit effect**; TMB, if anything, suppresses the association.

- **A2 — two-channel mediation (refutes the project's original novelty).** TMB-adjusted, 2,000
  bootstraps. APM channel mediates **28% [9,47]** (significant); the **cGAS-STING sensor channel is
  null (+1% [−3,6])**. The "sensor suppression is a new evasion channel" hypothesis — the reason the
  project felt exciting — **is dead**. The surviving channel (antigen presentation) **confirms**
  Foidart 2026 (WGD → MHC-I silencing); the contribution is causal quantification, not mechanism.

- **A2b — residual hunt (no clean new mechanism).** ~71% of WGD's effect is direct. The largest
  apparent mediators (Treg 37%, myeloid 34%) are **statistically circular** — co-linear
  immune-infiltration markers, not causal channels. The one mechanistically directional axis
  (stroma/TGF-β) does **not** mediate WGD (WGD tumours don't have more stroma). The residual reads as
  global immune desertification that bulk RNA cannot localise.

- **A5 — independent replication (THE HEADLINE FAILS).** Non-TCGA ICGC donors, same ShatterSeek
  exposure calls, independent PCAWG RNA. WGD effect: TCGA −0.339 (p=6.6e-7) → PCAWG −0.151 (p=0.17) →
  PCAWG solid-only −0.069 (p=0.56). All three benign explanations **tested and ruled out**:
  (1) tumour-type mix — restricting TCGA to the replication's histologies is *stronger* (−0.520,
  p=0.0017), i.e. same cancers, opposite answer; (2) lymphoid contamination — excluding it makes the
  estimate *smaller*; (3) data quality — positive controls (APM~cytolytic, sensor~cytolytic, purity)
  reproduce at TCGA strength. Power is adequate (SE~0.12 at n=260 would detect the TCGA magnitude at
  p<0.01). **The effect is cohort-dependent — most likely a TCGA-specific technical/selection or
  normalisation sensitivity — not established biology.**

---

## 3. What is and is not true, after everything

| Claim | Status |
|---|---|
| Chromothripsis→immune-cold is largely WGD acting through it (mutual adjustment) | **Real within TCGA**, TMB-robust; but see replication |
| The phenotype is a neoantigen-deficit (low-TMB) effect | **False** — WGD tumours have more mutations yet are colder |
| WGD evades immunity via cGAS-STING **sensor** suppression (the "new mechanism") | **Refuted** — sensor channel null in formal mediation |
| WGD → immune-cold runs through **antigen presentation** | Supported (28%), but **confirms** Foidart 2026 — not novel |
| A clean novel mediator explains the direct 71% | **No** — not stroma/EMT/WNT; apparent mediators are circular |
| WGD → immune-cold is established biology | **Not on this evidence** — fails independent replication |

---

## 4. Is anything publishable?

- **Yes, as a cautionary methods result — the honest angle.** A large fraction of the aneuploidy/
  WGD–immune literature is TCGA-only. Here is an effect that is robust *within* TCGA (p=6.6e-7,
  survives TMB, E-value 2.08) and yet **collapses in an independent cohort with every trivial
  explanation excluded**. That is a genuinely useful contribution: *how much do we trust
  single-cohort aneuploidy-immune associations?* It reframes the null as the finding.
- **The chromothripsis-vs-WGD dissociation** is worth stating regardless (the literature conflates
  them), but only *with* the replication caveat attached, not as a standalone positive claim.
- **The calcification label** (from the origin) is a small reusable resource/method (nobody had
  mined TCGA-Reports for calcification).
- **Not publishable:** as a positive WGD-mechanism paper, or as a new-mechanism (sensor) discovery.

---

## 5. Recurring methodological lessons (the real portfolio value)

1. **Replicate in an independent cohort before calling anything a flagship.** The strongest lesson:
   TMB-robust + E-value 2.08 *within TCGA* still failed to generalise. Statistical robustness inside
   one cohort is not replication. (This should be applied retroactively to the other flagships —
   ecDNA-immune, etc. — before they are pitched.)
2. **Bulk-RNA immune signals are weak and entangled** with proliferation, purity, infiltration; most
   "mediators" recovered this way are co-linear readouts, not mechanisms. Guard against the mediation
   circularity trap.
3. **Cheap empirical de-risking early saves weeks.** Three ideas were correctly killed for a few
   CPU-hours each; the survivor was found by anchoring to a PI's actual data, not by brainstorming.
4. **Honest negatives are the deliverable, not a failure** — provided they are reported as such and
   the salvageable byproducts (resources, cautionary results) are kept.

---

## 6. Recommendations

- **Do not submit as a positive finding.** Either shelve, or write the short **cautionary
  cohort-dependence** note (which is real and useful).
- **Retire the "sensor suppression" framing entirely** — it is refuted and the README/RESULTS
  already say so; keep it that way.
- **Apply the replication discipline to the other flagships** before any PI pitch.
- **Keep the collaboration hook honest:** the unresolved direct effect is a legitimate single-cell/
  spatial question for the group whose calls this is built on — but pitched as "here is a
  cohort-dependent bulk signal that needs single-cell resolution," not as an established mechanism.
- **Preserve** the calcification label resource and the report-mining machinery for reuse.

*Prepared as a standing analysis for this project. If the project is revisited, start from §4/§6.*
