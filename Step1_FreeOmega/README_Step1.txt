STEP 1: FREE OMEGAS ESTIMATION
===============================

Paper: "Who Bears the Brunt of Recessions?" (Bognanni, Hall, Kudlyak, Wiczer)
Step:  1 of 4 — Year-by-year estimation with free type shares

MODEL
-----
- 5 types: 3 mover (High-N, High-U, High-E) + 2 polar (All-N, All-E)
- 3 mover types: full 4x4 transition matrix estimated per year
- 2 polar types: hardwired (All-N = always OLF, All-E = always employed)
- Type shares (omegas) estimated freely each year (sum to 1)
- Data: CPS 1976-2024, prime-age men, 8-month activity sequences (6561 paths)
- Missing years: 1977, 1993 (no 8-month match possible)
- Objective: maximum likelihood (negative log-likelihood minimized)

CONFIG TOGGLES (in run_step1.m, overriding config_Paper2.m)
-----------------------------------------------------------
  cfg.mode              = 'free_omega'
  cfg.Theta             = 3           (3 mover types estimated)
  cfg.ThetaTot          = 5           (+ 2 polar = 5 total)
  cfg.omegaEstimated    = true        (omegas are free parameters)
  cfg.allE_estimated    = false       (All-E is hardwired, not estimated)
  cfg.useRegularization = false       (no regularization penalties)

PARAMETERS PER YEAR
-------------------
  36 transition matrix entries (3 types x 12 params each)
   5 type shares (omegas, constrained to sum to 1)
  = 41 free parameters

HOW TO RUN
----------
1. Open MATLAB
2. cd to the Paper2_Rewrite folder
3. Run:  run_step1
   This will:
   - Estimate all 47 years in parallel (parfor)
   - Export results to Excel (Step1_FreeOmega_Results.xlsx)
   - Generate 6 publication-quality figures (.png + .fig)
   Estimated time: ~30 minutes on 16-core laptop

OUTPUT FILES
------------
  Step1_FreeOmega/
    step1_1.mat ... step1_47.mat    — per-year result structs
    convergence_step1.csv           — convergence diagnostics
    convergence_step1.mat           — same, in .mat format
    Step1_FreeOmega_Results.xlsx    — Excel with 5 sheets:
        Omegas      : type shares over time
        Ergodic     : ergodic distribution by type
        Rates       : U-rate, OLF-rate, E-rate by type + aggregate
        Parameters  : all 12 transition params per type
        Convergence : exit flags, objective values
    Fig1_omegas.png                 — all type shares (cf. slide 51)
    Fig2_omegas_individual.png      — High-U, High-E, All-E panels (cf. slide 52)
    Fig3_urate_mover.png            — U-rate by mover type (cf. slide 53)
    Fig4_ergodic_by_type.png        — 4-panel ergodic distributions
    Fig5_urate_all.png              — U-rate by type + aggregate
    Fig6_olf_all.png                — OLF rate by type + aggregate

TYPE LABELING (post-estimation)
-------------------------------
Types are labeled after estimation based on ergodic distribution:
  1. Sort by ergodic E mass (states 3+4): lowest = High-N, highest = High-E
  2. Sub-sort the bottom two by U mass: higher U mass = High-U
  3. States ordered so p33 <= p44 (short-term E vs long-term E)

KEY RESULTS TO CHECK
--------------------
- High-U unemployment rate should be procyclical (slide 53)
- Type shares should be roughly: All-E ~0.65-0.80, rest smaller
- Exit flag should be 1 for most years
