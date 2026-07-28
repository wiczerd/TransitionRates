# Step 1: Free Omega Estimation — Working Log

## 2026-03-11: Initial Diagnostic Review

### Run Summary
- **Code used:** Original MATLAB code (`ModelFitFuncTimeSeries_NewTrMBoot_20250919_V2.m`) via wrapper `run_step1_original.m`
- **Data:** CPS 1976–2024, prime-age men, 8-month activity sequences (6561 paths, 47 years)
- **Model:** 5 types (3 mover + 2 polar), 41 free parameters per year (36 transition + 5 omegas)
- **Results location:** `Paper2_Rewrite/Step1_FreeOmega/`
- **Note:** The folder contains a `run_step1_log.txt` (from an earlier run at 16:27) and `convergence_step1.csv` (from the later, better run at 19:15). The .mat files and figures correspond to the later run (the CSV).

### Optimizer Structure (from `run_step1_original.m`)
The code runs up to 3 optimizer passes per year:
1. **Pass 1:** `fmincon` interior-point, forward finite differences
2. **Pass 2** (if Pass 1 gave exitflag=2 or firstorderopt > 0.001): `fmincon` SQP, central differences, L-BFGS Hessian
3. **Pass 3** (if Pass 2 still gave exitflag=2 or 0): `fmincon` interior-point, central differences

### Convergence
All 47/47 years reached exitflag=1 in the final run. However, **6 years have poor first-order optimality**, meaning even after all 3 passes the optimizer didn't converge below the 0.001 tolerance:

| Year | FirstOrderOpt | Fval   | Comment                        |
|------|---------------|--------|--------------------------------|
| 1984 | 0.358         | 1.7229 | Should be < 0.001              |
| 1987 | 0.448         | 1.5223 | Should be < 0.001              |
| 2003 | **0.751**     | 1.9214 | Worst — far from a local min   |
| 2010 | 0.347         | 2.1151 |                                |
| 2011 | 0.475         | 2.0292 |                                |
| 2013 | 0.523         | 1.8926 |                                |

The remaining 41 years all have FirstOrderOpt < 0.001.

### Key Finding: Estimates Are Too Noisy (Non-Convexity Problem)

The estimated parameter series fluctuate far more than the underlying data warrants. Since CPS data changes smoothly year-to-year, wild jumps in estimated parameters indicate the optimizer is landing in **different local optima** in different years. Evidence:

- **Type shares (Fig1):** All-E is roughly stable (~55–75%), but mover type shares (High-E, High-U, High-N) jump wildly between adjacent years.
- **Unemployment rates (Fig3):** High-U bounces between 20–80% with no coherent cyclical pattern visible through the noise. High-N similarly erratic.
- **Ergodic distributions (Fig4):** All four panels (OLF, U, E-short, E-long) are extremely noisy. Year-to-year variation in ergodic shares is far larger than any plausible economic change.

### Suspected Type-Label Switching

The noisy spikes in Fig3 and Fig4 suggest that in some years the "High-U" and "High-N" types may be swapping labels. The post-estimation labeling algorithm (sort by ergodic E mass, sub-sort by U mass) can flip the ordering when applied to slightly different optima — creating the appearance of wild jumps even if the underlying types are actually similar.

### Assessment

The run finished and gives a first look at the results, but **the estimates are not publication-quality** due to:
1. Multiple local optima across years (noisy parameter paths)
2. Six years not properly converged
3. Possible type-label switching

### Strategies to Explore (not yet attempted)

- **Warm-starting** each year from its neighbor's solution (keep optimizer in same basin)
- **Multi-start** with many more random initial guesses per year
- **Continuity penalty** across adjacent years (smooth regularization)
- **Fixing the 6 poorly-converged years** first, then reassessing

---

### Clarifications Made During Discussion

**Log file vs CSV discrepancy:** Initially described as a "two-stage process" within one run. After checking timestamps, we confirmed these are from **two separate runs** (~3 hours apart). The `run_step1_log.txt` is from an earlier run (16:27, 42/47 exitflag=1). The `convergence_step1.csv`, all .mat files, and all figures are from the later, better run (19:15, 47/47 exitflag=1). The earlier log was never overwritten.

**"SQP refinement" terminology:** This is accurate — the code at lines 154–177 of `run_step1_original.m` implements a conditional multi-pass strategy. Pass 2 uses the SQP algorithm with central differences and L-BFGS Hessian. But the two different result sets (log vs CSV) are from two different runs of this same multi-pass code, not from different passes within one run.

### Session Logistics

- **Memory files** stay in `.claude` (auto-loaded each session), not moved to Dropbox
- **This working log** (`step1_freeomega_log.md`) lives with the results in `Step1_FreeOmega/` — one log per phase, not daily
- Claude Code runs locally; no background work happens when the computer is off
- Next session: pick up from the strategies listed above

---

## 2026-03-12: Multi-Start Script + Robustness Strategy

### Discussion of Cyclical Patterns

Marianna reviewed the Fig1 and Fig3 results and **disagrees that the estimates are purely noise**. Her reading:

- All-E and All-N shares are smooth — as expected
- Mover type shares bounce but in an **economically logical** direction: around recessions, High-U share increases while High-E share declines
- The unemployment rate of High-U **declines** in recessions — this makes sense because if omegas are free, the model re-allocates more individuals into the High-U group in recessions but dilutes their unemployment rate
- This cyclical pattern is the **foundation for Phases 2 and 3** of the paper: type shares should be smooth (fixed omegas), and transition matrices should absorb the time-series variation

### Type Labeling

Current criterion: sort by ergodic E mass, sub-sort by U mass. Marianna confirms this is the right economic criterion — she does not want to call "High-E" a type that doesn't have the highest E share in ergodic among movers. Label-switching concern is deferred until after improving the optimization. May revisit if improved estimates still show erratic jumps.

### Decision: Fix All 47 Years with Multi-Start

Rather than fixing only the 6 problematic years, we decided to **re-run all 47 years** with multi-start to potentially improve solutions everywhere.

### Script Created: `rerun_all_multistart.m`

**Location:** `Paper2_Rewrite/rerun_all_multistart.m`
**Output:** `Paper2_Rewrite/Step1_FreeOmega_Improving/` (new subfolder; old results in `Step1_FreeOmega/` preserved untouched)

**Two-phase approach per year:**

**Phase A — Quick screening** (MaxIter=3000, OptTols=0.01): Run 7 candidates to find the best basin:
1. Old Step1 solution (self warm-start)
2. Neighbor t−1 solution
3. Neighbor t+1 solution
4–5. Jittered old solution (Gaussian noise σ=0.03, omegas renormalized)
6–7. Random initial guesses from seeds 200, 300 (original used seed 100)

**Phase B — Full refinement**: Same 3-pass optimizer as original (interior-point → SQP → interior-point central) on the best candidate from Phase A.

**Everything else identical** to `run_step1_original.m`: objective function, constraints, post-processing, output format.

**Additional output:** `convergence_multistart.csv` with old vs. new fval comparison for every year.

**Estimated runtime:** ~2–2.5 hours with 6 workers.

**Status:** Script written, NOT YET RUN.

### Broader Robustness Strategy for the Paper

Discussed how to convince a referee that the solution is the global optimum. Strategies in priority order:

1. **Monte Carlo validation** (essential) — Simulate data from estimated parameters, re-estimate, check parameter recovery. Proves the procedure works when the truth is known.

2. **EM algorithm** (highly recommended) — The classical method for mixture models (Dempster, Laird & Rubin 1977). E-step: posterior type probabilities given data. M-step: update transition matrices and omegas. If EM and fmincon agree, that's very strong evidence. EM also has monotonic likelihood improvement and is less sensitive to starting values.

3. **Multi-start fmincon** (already doing) — Report number of starting points, show best solution is robust across starts.

4. **BIC for number of types** — Justify ThetaTot=5 vs alternatives (4 or 6). Quick to compute, addresses obvious referee question.

5. **Profile likelihood** — For key parameters (e.g., omega_HighU), fix at a grid of values, optimize over rest, plot. Shows whether optimum is sharp or flat. A few representative years suffice.

6. **Bayesian MCMC** (gold standard but high effort) — Full posterior sampling, handles multi-modality, gives uncertainty quantification.

### Next Steps

1. Run `rerun_all_multistart.m` and compare results with old Step1
2. Generate figures for the improved estimates
3. Begin implementing Monte Carlo validation
4. Begin implementing EM algorithm
5. Reassess type-labeling and smoothness after seeing improved results

---

## 2026-03-13: Multi-Start Run + EM Discussion

### Decision: Multi-Start First, Then EM

Discussed whether to proceed with (1) multi-start for each year or (2) EM algorithm. Decision: **run multi-start first** because:
- Script `rerun_all_multistart.m` already written and ready
- Provides immediate improvement baseline with zero new coding
- Results inform whether EM is needed for confirmation or improvement

### Multi-Start Run — COMPLETED

**Status:** Completed successfully
**Output folder:** `Paper2_Rewrite/Step1_FreeOmega_Improving/`
**Workers:** 6 parallel
**Total runtime:** 119.4 minutes (~2 hours)
**Results file:** `Step1_FreeOmega_Improving/convergence_multistart.csv`

### Final Results: 47/47 exitflag=1, 12/47 improved

All 47 years converged (exitflag=1). 12 years found a better solution than the original run.

#### The 6 Previously Problematic Years — 5/6 Substantially Improved

| Year | Old fval | New fval | Delta | Old 1stOrdOpt | New 1stOrdOpt | Best Candidate |
|------|----------|----------|-------|---------------|---------------|----------------|
| **1984** | 1.7229 | 1.7020 | **-0.0209** | 0.358 | 0.00086 | 2 (neighbor t-1) |
| **1987** | 1.5223 | 1.5097 | **-0.0126** | 0.448 | 0.001 | 3 (neighbor t+1) |
| **2003** | 1.9214 | 1.9225 | +0.0011 | 0.751 | 0.001 | 6 (random) |
| **2010** | 2.1151 | 2.0968 | **-0.0183** | 0.347 | 0.001 | 2 (neighbor t-1) |
| **2011** | 2.0292 | 2.0141 | **-0.0151** | 0.475 | 0.00088 | 5 (jittered) |
| **2013** | 1.8926 | 1.8806 | **-0.0119** | 0.523 | 0.001 | 3 (neighbor t+1) |

All 6 now have firstorderopt ≤ 0.001 (previously 0.35–0.75). 2003 didn't improve fval but convergence is now proper.

#### Other Notable Improvements

| Year | Old fval | New fval | Delta | Best Candidate |
|------|----------|----------|-------|----------------|
| 1989 | 1.5192 | 1.5114 | -0.0078 | 3 (neighbor t+1) |
| 1992 | 1.7274 | 1.7214 | -0.0061 | 2 (neighbor t-1) |
| 2009 | 2.1366 | 2.1312 | -0.0054 | 1 (old, refined better) |
| 2015 | 1.7685 | 1.7633 | -0.0052 | 2 (neighbor t-1) |
| 2018 | 1.6520 | 1.6473 | -0.0047 | 2 (neighbor t-1) |
| 2024 | 1.5785 | 1.5763 | -0.0022 | 2 (neighbor t-1) |

#### New Problem: Year 2024

Year 2024 has **firstorderopt = 0.712** — it did NOT converge well despite improving fval slightly. This is a new problem (old run had it converged). Needs targeted re-run.

#### Key Takeaway

Neighbor solutions (candidates 2 and 3) were the most valuable starting points — they won for most improved years. This confirms year-to-year solutions share basins, and warm-starting from neighbors navigates the non-convex landscape effectively.

### EM Algorithm Explained

Discussed what EM is and why it's relevant:
- **E-step:** Given current params, compute posterior probability each path belongs to each type (Bayes' rule)
- **M-step:** Given posteriors, update omegas (weighted average) and transition matrices (smaller 12-param optimizations per type)
- **Advantages:** Monotonic likelihood improvement, simpler sub-problems, less sensitive to starting values, natural for mixture models
- **Value:** If EM and fmincon agree → strong evidence of global optimum (convincing for referees)
- **Limitation:** Can also get stuck in local optima; M-step for transition matrices may not be closed-form

### Next Steps (for next session)

1. Fix Year 2024 convergence issue
2. Generate figures for improved estimates — compare smoothness with old results
3. Assess whether improvements change the cyclical patterns
4. Decide on EM implementation based on visual inspection of improved results

---

## 2026-03-14: Year 2024 Fix, EM Algorithm, Ergodic Discussion

### Year 2024 Fix — DONE
- Script: `fix_year2024.m` — 20 candidates, lb=1e-6 (not 0), top 3 refined
- **fval: 1.5785 → 1.5644** (improvement of 0.014, largest single-year gain)
- **firstorderopt: 0.712 → 0.00008** (excellent convergence)
- Winner: Neighbor 2019 (not adjacent years 2022/2023)
- Solution shifted dramatically: High-E=0.63/All-E=0.05 (was 0.28/0.58)
- No parameters at bounds (lb=1e-6 fix worked)
- Saved to `Step1_FreeOmega_Best/step1_47.mat` and `Step1_FreeOmega_Improving/step1_47.mat`

### Best-of-Both Figures — DONE
- Script: `build_best_of_both.m` + `figures_compare_multistart.m`
- Results in `Step1_FreeOmega_Best/Figures/` and `Step1_FreeOmega_Improving/Figures/`
- 13/47 years use multi-start, 34/47 keep original
- Multi-start improved worst years but fundamental noisiness of mover types persists
- Supports moving to fixed omegas (Phases 2/3)

### EM Algorithm — COMPLETED
- Script: `run_em_step1.m` — EM with CPS gap-aware M-step, 5 starts/year
- Results in `Step1_FreeOmega_EM/` — runtime 50 min, 6 workers
- **fmincon wins 37/47 years**; EM wins 10/47 with tiny improvements
- Key finding: EM eliminates All-N polar type in many years (omega→0)
- 21/47 years methods agree (|delta fval|<0.001); 26 disagree
- Low omega correlations between methods (0.20-0.52)
- EM underperforms likely due to approximate M-step (ergodic constraint)
- Figures in `Step1_FreeOmega_EM/Figures/`

### Ergodic Initial Distribution Discussion
- In Paper 1 (cross-section): assumption reasonable
- In Paper 2 Phase 1 (free omegas): tension exists but omegas compensate
- In Paper 2 Phases 2/3 (fixed omegas): tension worse — could dampen cyclicality
- EM finding (killing All-N) related to initial distribution sensitivity
- **Alternatives discussed:** free initial dist (+9 params), previous year's ergodic, weighted avg, keep + acknowledge
- **Decision: implement free initial distribution as robustness check**

### Next Steps
1. Implement free initial distribution robustness check
2. Monte Carlo validation
3. BIC for number of types
4. Move to Phases 2/3
