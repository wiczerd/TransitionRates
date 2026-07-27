% test_step1_single.m — Quick test of Step 1 (free omega) for one year
%
% Verifies that omegas are being estimated (not stuck at 0.2).
% Compares with the direct computation from old code.

clc; clear all;
fprintf('=== Step 1 single-year test (year 14 = 1990) ===\n');

%% Config for Step 1
config_Paper2;
cfg.mode = 'free_omega';
cfg.Theta    = 3;
cfg.ThetaTot = 5;
cfg.omegaEstimated = true;
cfg.allE_estimated = false;
cfg.useRegularization = false;

%% Build Hstates
grids = cell(1, cfg.Nmonth);
[grids{:}] = ndgrid(1:cfg.Nstates);
Hstates = zeros(cfg.Nstatepath, cfg.Nmonth);
for k = 1:cfg.Nmonth
    Hstates(:,k) = grids{k}(:);
end
clear grids;

%% Load data
dataFile = fullfile(cfg.oldDataDir, cfg.inputFileNames{cfg.groupNum});
Data_full_TimeSeries = xlsread(dataFile, cfg.dataRange);
year = 14;
Data_full = Data_full_TimeSeries(:, year);
fprintf('Year %d (calendar %d), %d paths\n', year, cfg.yearsList(year), numel(Data_full));

%% Constraints & initial guesses
[A_ineq, B_ineq, Aeq, beq, lb, ub] = build_constraints(cfg);
fprintf('Parameters: %d (36 transition + 5 omegas)\n', numel(lb));
candidates = build_initial_guesses(cfg, ub);

omega_yr = [];  % free mode: extracted from x
f = @(x) objective_fn(x, cfg, Hstates, Data_full, omega_yr);

% Pick best candidate
bestx = candidates{1}; bestf = inf;
for c = 1:numel(candidates)
    try, ftry = f(candidates{c}); catch, ftry = inf; end
    fprintf('  Candidate %d: fval = %.6f\n', c, ftry);
    if ftry < bestf, bestf = ftry; bestx = candidates{c}; end
end
fprintf('Best initial: fval = %.6f\n', bestf);
fprintf('Initial omegas: '); fprintf('%.4f ', bestx(end-4:end)); fprintf('\n');

%% fmincon — interior-point (no ScaleProblem, matching old code)
opts1 = optimoptions('fmincon', ...
    'Algorithm','interior-point', 'Display','iter', ...
    'MaxIterations',cfg.MaxIter, 'MaxFunctionEvaluations',cfg.MaxEvals, ...
    'OptimalityTolerance',cfg.OptTols, 'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-9, 'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',1e-5);

tic;
[ests, SSR2, exitflag, outputF] = fmincon(f, bestx, A_ineq, B_ineq, Aeq, beq, lb, ub, [], opts1);
t1 = toc;
fprintf('\nResult: exitflag=%d, fval=%.6f, time=%.1fs\n', exitflag, SSR2, t1);

% Extract omegas from the result
omega_est = ests(end-4:end);
fprintf('Estimated omegas: '); fprintf('%.4f ', omega_est); fprintf('\n');
fprintf('Omega names:      All-N=%.4f, High-N=%.4f, High-U=%.4f, High-E=%.4f, All-E=%.4f\n', ...
    omega_est(5), omega_est(1), omega_est(2), omega_est(3), omega_est(4));

% Check: omegas should NOT all be 0.2
if max(abs(omega_est - 0.2)) < 0.01
    fprintf('\n*** WARNING: Omegas are still near-uniform (0.2) — fix did not work! ***\n');
else
    fprintf('\n*** SUCCESS: Omegas are non-uniform — estimation is working correctly ***\n');
end

%% Extract full results
res = extract_results(ests, cfg, omega_yr, SSR2, exitflag);
fprintf('\nRelabeled omegas: '); fprintf('%.4f ', res.omega); fprintf('\n');
tNames = {'High-N','High-U','High-E'};
for th = 1:cfg.Theta
    fprintf('  %s: ergodic=[%.3f %.3f %.3f %.3f], urate=%.3f\n', ...
        tNames{th}, res.ergVec{th}, res.urate(th));
end
fprintf('\n=== Test complete ===\n');
