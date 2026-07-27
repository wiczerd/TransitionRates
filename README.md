# Paper 2 Rewrite — MATLAB Code

All scripts created by Claude Code to support the estimation and analysis in
"Who Bears the Brunt of Recessions?" (Bognanni, Hall, Kudlyak, Wiczer).

Working directory: `2026 Matlab code/Paper2_Rewrite/`

---

## Shared Utilities

| File | Description |
|------|-------------|
| `Ergodic.m` | Computes the ergodic (stationary) distribution of a transition matrix. |
| `build_constraints.m` | Builds inequality constraints (A, b) and bounds (lb, ub) for fmincon: row-sum constraints on TM parameters. |
| `build_initial_guesses.m` | Generates multiple starting points for multi-start optimization from neighbors, random draws, and perturbations. |
| `config_Paper2.m` | Central configuration: data paths, year list, group number, Theta/ThetaTot, Hstates matrix. Sourced by all estimation scripts. |
| `extract_results.m` | Extracts ergodic distributions, unemployment rates, and omegas from saved `.mat` result files. |

---

## Phase 1: Free Omegas (Year-by-Year Estimation)

### Objective Functions

| File | Description |
|------|-------------|
| `objective_fn.m` | Phase 1 objective function (negative log-likelihood) for the mixture model with free omegas and ergodic initial distribution. Wraps the original `ModelFitFuncTimeSeries` logic. |
| `ModelFitFuncTimeSeries_FreeInitDist.m` | Phase 1 objective with **free initial distribution** (50 params: 36 TM + 5 omegas + 9 init dist). Used by `run_step1_freeinit.m` and `run_em_step1_freeinit.m`. |

### Estimation Scripts

| File | Description |
|------|-------------|
| `run_step1.m` | Phase 1 main estimation: multi-start fmincon for all 47 years. 4-pass optimization (interior-point, SQP, central differences, polish). Results → `Step1_FreeOmega/`. |
| `run_step1_original.m` | Runs the user's original Phase 1 code on the rewritten data pipeline for comparison/replication. Results → `Step1_FreeOmega/`. |
| `run_estimation.m` | Earlier version of Phase 1 estimation (superseded by `run_step1.m`). |
| `rerun_all_multistart.m` | Multi-start re-estimation for all 47 years with 9 candidates per year (neighbor warm-starts + random + perturbations). Results → `Step1_FreeOmega_Improving/`. |
| `run_step1_freeinit.m` | Phase 1 with free initial distribution using fmincon multi-start. 50 params per year. Results → `Step1_FreeOmega_FreeInitDist/`. |
| `run_em_step1.m` | Phase 1 EM algorithm with ergodic initial distribution. Gap-aware M-step handles the 8-month CPS panel gap. Results → `Step1_FreeOmega_EM/`. |
| `run_em_step1_freeinit.m` | Phase 1 EM algorithm with **free initial distribution**. Exact closed-form M-step for pi (no approximation needed). Best Phase 1 results. Results → `Step1_FreeOmega_FreeInitDist_EM/`. |

### Fix Scripts

| File | Description |
|------|-------------|
| `fix_year2024.m` | Targeted fix for year 2024 (poorly converged): 20 candidates with lb=1e-6. |
| `fix_year2024_v2.m` | Second fix for year 2024 type swap: 30 candidates, non-swapped solution. Saved to `Step1_FreeOmega_Best/`. |
| `fix_freeinit_worse_years.m` | Fixes 8 years where free-init fmincon was worse than ergodic baseline, using targeted multi-start. |
| `build_best_of_both.m` | Combines original + multi-start results: picks min(fval) per year. Results → `Step1_FreeOmega_Best/`. |

### Figures & Diagnostics

| File | Description |
|------|-------------|
| `figures_step1.m` | Generates Phase 1 summary figures (omegas, U-rates, ergodic distributions) from `Step1_FreeOmega/`. |
| `figures_step1_orig.m` | Same figures from the original code replication results. |
| `figures_compare_multistart.m` | 7 comparison figures: original vs multi-start vs best-of-both. Saved to `Step1_FreeOmega_Best/Figures/`. |
| `figures_compare_em.m` | 5 comparison figures: fmincon vs EM (both ergodic). Saved to `Step1_FreeOmega_EM/Figures/`. |
| `figures_compare_freeinit.m` | Comparison figures: fmincon ergodic vs fmincon free-init. Saved to `Step1_FreeOmega_FreeInitDist/Figures/`. |
| `figures_compare_em_freeinit.m` | 7 comparison figures across 3 methods (fmincon ergodic, fmincon free-init, EM free-init). Saved to `Step1_FreeOmega_FreeInitDist_EM/Figures/`. |
| `check_type_labeling.m` | Compares E-U-N vs By-E labeling criteria across 47 Phase 1 years. |

### Testing

| File | Description |
|------|-------------|
| `test_step1_single.m` | Test script: runs Phase 1 estimation for a single year to verify setup. |
| `test_year14.m` | Test script: compares our Phase 1 year-14 result against `guess_14.mat` warm-start. |
| `compare_year14.m` | Detailed comparison of year-14 solution vs original code. |
| `debug_gradient.m` | Finite-difference gradient check for the Phase 1 objective function. |
| `export_step1.m` | Exports Phase 1 results to CSV/Excel for external analysis. |
| `export_step1_orig.m` | Exports original-code Phase 1 results. |
| `plot_results.m` | General-purpose plotting of Phase 1 results. |

---

## Phase 2: Omega Trends

| File | Description |
|------|-------------|
| `compute_omega_trends.m` | Computes 5 trend candidates for type shares (2 HP filters + 3 linear trends), all normalized to sum=1. Saves to `Phase2_OmegaTrends/`. Generates per-type comparison figures. |
| `create_omega_fixed.m` | Creates `omega_fixed.mat` from Phase 2 trend 3 (linear excl 2020-2021) for use in Phase 3. |

---

## Phase 3: Fixed Omegas, Time-Varying Transition Matrices

### Objective Functions

| File | Description |
|------|-------------|
| `ModelFitFunc_FixedOmega_AllE.m` | Phase 3 objective: fixed omegas, All-E as estimated type (Theta=4). Includes regularization (lam=2.0, tau=0.05, m=0.02) to prevent type swapping. 48 free params. Ergodic initial distribution. |
| `ModelFitFunc_FixedOmega_AllE_NoReg.m` | Same as above but without regularization. Used for testing. |
| `ModelFitFunc_FixedOmega_AllE_FreeInit.m` | Phase 3 objective with **free initial distribution**. 60 params (48 TM + 12 init dist). Regularization based on ergodic distributions (not init dist). |

### Estimation Scripts

| File | Description |
|------|-------------|
| `run_step3_fixedomega.m` | Phase 3 main estimation: 4-pass optimization with parfor parallelism. Loads `omega_fixed.mat`. Results → `Step3_FixedOmega/`. |
| `run_step3_sequential.m` | Phase 3 sequential warm-starting: forward pass (1→47), backward pass (47→1), best-of-three. Reduces noise from basin-switching. Results → `Step3_FixedOmega_Sequential/`. |
| `run_step3_em_freeinit.m` | Phase 3 sequential EM with free initial distribution. Forward+backward+best-of-three. 3 EM starts per year + fmincon refinement. Results → `Step3_FixedOmega_FreeInit_EM/`. |

### Fix & Diagnostic Scripts

| File | Description |
|------|-------------|
| `fix_step3_nonconverged.m` | Fixes 4 non-converged Phase 3 years (1987, 2022, 2023, 2024) with 30 candidates each using Phase 3 neighbor warm-starts. |
| `check_step3_labeling.m` | Diagnoses Phase 3 type labeling: compares By-E vs E-U-N criteria. |
| `fix_step3_labeling.m` | Re-sorts 14 Phase 3 years from By-E to E-U-N labeling criterion. |
| `test_step3_year20.m` | Test script: runs Phase 3 for a single year to verify setup. |

### Figures

| File | Description |
|------|-------------|
| `figures_step3.m` | 7 Phase 3 figures (omegas, E/U/N mass, U-rates, stacked bars). Saved to `Step3_FixedOmega/Figures/`. |
| `figures_step3_sequential.m` | 7 comparison figures: sequential vs existing Phase 3. Saved to `Step3_FixedOmega_Sequential/Figures/`. |
| `figures_step3_paper.m` | **13 clean paper-ready figures** (sequential results only, with transition rates). `paper_fig1` through `paper_fig13`. Saved to `Step3_FixedOmega_Sequential/Figures/`. |
| `decomposition_table.m` | Computes recession decomposition table + 2 figures (stacked bar, peak/trough). Saves CSV + figures to `Step3_FixedOmega_Sequential/`. |

---

## Post-Processing

| File | Description |
|------|-------------|
| `compute_fit_stats.m` | Computes model fit statistics (data entropy, KL divergence, pseudo-R²) from Phase 3 sequential results. Saves to `Step3_FixedOmega_Sequential/fit_stats.mat`. **Queued to run after overnight EM job.** |

---

## Output Directories

| Directory | Contents |
|-----------|----------|
| `Step1_FreeOmega/` | Phase 1 results (original run, 6 years poorly converged) |
| `Step1_FreeOmega_Improving/` | Phase 1 multi-start results (12/47 improved) |
| `Step1_FreeOmega_Best/` | Best-of-both: min(fval) per year from original + multi-start |
| `Step1_FreeOmega_EM/` | Phase 1 EM algorithm results |
| `Step1_FreeOmega_FreeInitDist/` | Phase 1 fmincon with free initial distribution |
| `Step1_FreeOmega_FreeInitDist_EM/` | Phase 1 EM with free initial distribution (**best Phase 1**) |
| `Phase2_OmegaTrends/` | Trend estimates and comparison figures |
| `Step3_FixedOmega/` | Phase 3 initial results (parfor, multi-start) |
| `Step3_FixedOmega_Sequential/` | Phase 3 sequential warm-starting results (**main Phase 3**) |
| `Step3_FixedOmega_FreeInit_EM/` | Phase 3 sequential EM with free init dist (running) |

---

## Data Files (not created by Claude, used as inputs)

| File | Description |
|------|-------------|
| `DataTimeSeriesPM_1976_2023.xlsx` | Activity path frequencies, prime-age men, columns J-BD (47 years) |
| `DataTimeSeriesYW.xlsx`, `DataTimeSeriesPW.xlsx`, `DataTimeSeriesYM.xlsx` | Other demographic groups |
| `omega_fixed.mat` | Phase 2 trend type shares (created by `create_omega_fixed.m`) |
| `guess_14.mat` | Warm-start solution from year 14 (original code) |

---

## Execution Order

A typical full replication proceeds as:

1. `run_step1.m` → Phase 1 baseline
2. `rerun_all_multistart.m` → Phase 1 multi-start
3. `build_best_of_both.m` → combine best results
4. `fix_year2024_v2.m` → fix 2024 type swap
5. `run_em_step1.m` → Phase 1 EM (robustness)
6. `run_step1_freeinit.m` → Phase 1 free init dist (robustness)
7. `run_em_step1_freeinit.m` → Phase 1 EM + free init (robustness)
8. `compute_omega_trends.m` → Phase 2 trends
9. `create_omega_fixed.m` → create omega_fixed.mat
10. `run_step3_fixedomega.m` → Phase 3 baseline
11. `fix_step3_nonconverged.m` → fix non-converged years
12. `fix_step3_labeling.m` → fix type labeling
13. `run_step3_sequential.m` → Phase 3 sequential (main results)
14. `run_step3_em_freeinit.m` → Phase 3 EM + free init (robustness)
15. `figures_step3_paper.m` → paper figures
16. `decomposition_table.m` → recession decomposition
17. `compute_fit_stats.m` → model fit statistics
