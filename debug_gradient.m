% debug_gradient.m — Manually check gradient w.r.t. omegas at initial point
%
% Verifies whether the objective function actually changes when omegas move.

clc; clear all;
fprintf('=== Debug: gradient check at initial point ===\n\n');

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

%% Build initial guess (same as build_initial_guesses candidate 1)
[~, ~, ~, ~, ~, ub] = build_constraints(cfg);
candidates = build_initial_guesses(cfg, ub);
x0 = candidates{1};  % random start with uniform omegas

omega_yr = [];
f = @(x) objective_fn(x, cfg, Hstates, Data_full, omega_yr);

f0 = f(x0);
fprintf('Base objective at x0: %.10f\n', f0);
fprintf('Initial omegas: '); fprintf('%.4f ', x0(end-4:end)); fprintf('\n\n');

%% Manual perturbation of omegas
% Shift omega: increase All-E, decrease High-N (maintaining sum=1)
h_vals = [0.01, 0.05, 0.10, 0.20];
fprintf('--- Perturbation: increase omega_AllE, decrease omega_HighN ---\n');
fprintf('%10s %15s %15s\n', 'delta', 'f(x_perturbed)', 'f-f0');
for i = 1:numel(h_vals)
    h = h_vals(i);
    x1 = x0;
    x1(end-1) = x1(end-1) + h;  % omega_AllE (index 4 from end = end-1)
    x1(end-4) = x1(end-4) - h;  % omega_HighN (index 1 from end = end-4)
    f1 = f(x1);
    fprintf('%10.4f %15.10f %15.10f\n', h, f1, f1-f0);
end

fprintf('\n--- Perturbation: increase omega_AllE, decrease others equally ---\n');
for i = 1:numel(h_vals)
    h = h_vals(i);
    x1 = x0;
    x1(end-1) = x1(end-1) + h;         % omega_AllE += h
    x1(end-4) = x1(end-4) - h/4;       % others -= h/4 each
    x1(end-3) = x1(end-3) - h/4;
    x1(end-2) = x1(end-2) - h/4;
    x1(end)   = x1(end)   - h/4;
    f1 = f(x1);
    fprintf('%10.4f %15.10f %15.10f\n', h, f1, f1-f0);
end

%% Numerical gradient w.r.t. each omega (one-sided, h=1e-5)
fprintf('\n--- Finite-difference gradient w.r.t. omegas (forward, h=1e-5) ---\n');
h = 1e-5;
for k = 0:4
    idx = numel(x0) - 4 + k;  % omega indices: end-4 to end
    x_plus = x0;
    x_plus(idx) = x_plus(idx) + h;
    f_plus = f(x_plus);
    grad_k = (f_plus - f0) / h;
    fprintf('  omega(%d) [idx=%d]: df/domega = %.6f\n', k+1, idx, grad_k);
end

%% Check: what if we start with non-uniform omegas?
fprintf('\n--- Testing with non-uniform initial omegas ---\n');
x_nonunif = x0;
% Set omegas to roughly what the paper expects
x_nonunif(end-4) = 0.05;   % mover1 (will be relabeled)
x_nonunif(end-3) = 0.05;   % mover2
x_nonunif(end-2) = 0.10;   % mover3
x_nonunif(end-1) = 0.70;   % All-E
x_nonunif(end)   = 0.10;   % All-N
f_nonunif = f(x_nonunif);
fprintf('Objective with non-uniform omegas: %.10f (vs uniform: %.10f)\n', f_nonunif, f0);
fprintf('Difference: %.10f\n', f_nonunif - f0);

%% Try running fmincon from non-uniform start
fprintf('\n--- Running fmincon from non-uniform omegas ---\n');
[A_ineq, B_ineq, Aeq, beq, lb, ub] = build_constraints(cfg);
opts = optimoptions('fmincon', ...
    'Algorithm','interior-point', 'Display','iter', ...
    'MaxIterations',cfg.MaxIter, 'MaxFunctionEvaluations',cfg.MaxEvals, ...
    'OptimalityTolerance',cfg.OptTols, 'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-9, 'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',1e-5);

[ests, SSR2, exitflag, outputF] = fmincon(f, x_nonunif, A_ineq, B_ineq, Aeq, beq, lb, ub, [], opts);
fprintf('Result: exitflag=%d, fval=%.6f, iters=%d\n', exitflag, SSR2, outputF.iterations);
fprintf('Estimated omegas: '); fprintf('%.4f ', ests(end-4:end)); fprintf('\n');

fprintf('\n=== Debug complete ===\n');
