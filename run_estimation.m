% run_estimation.m — Main estimation driver for Paper 2
%
% Usage:
%   1. Edit config_Paper2.m to set mode ('free_omega' or 'fixed_omega')
%      and sub-mode ('both_polar_hardwired' or 'allE_estimated')
%   2. Run this script

clc; clear all;

%% Load configuration
config_Paper2;
fprintf('=== Paper 2 Estimation: mode=%s', cfg.mode);
if ~cfg.omegaEstimated
    fprintf(', submode=%s', cfg.fixedOmegaSubmode);
end
fprintf(' ===\n');
fprintf('Theta=%d estimated types, ThetaTot=%d total types, %d years\n', ...
    cfg.Theta, cfg.ThetaTot, cfg.numYears);

%% Build Hstates matrix (vectorized with ndgrid, replacing 8-deep nested loops)
grids = cell(1, cfg.Nmonth);
[grids{:}] = ndgrid(1:cfg.Nstates);
Hstates = zeros(cfg.Nstatepath, cfg.Nmonth);
for k = 1:cfg.Nmonth
    Hstates(:,k) = grids{k}(:);
end
clear grids;

%% Load data (once for all years)
dataFile = fullfile(cfg.oldDataDir, cfg.inputFileNames{cfg.groupNum});
fprintf('Loading data from %s ...\n', dataFile);
Data_full_TimeSeries = xlsread(dataFile, cfg.dataRange);
fprintf('  Data size: %d paths x %d years\n', size(Data_full_TimeSeries));

%% Load fixed omegas if needed
omega_fixed_all = [];
if ~cfg.omegaEstimated
    fprintf('Loading fixed omegas from %s ...\n', cfg.omegaFile);
    tmp = load(cfg.omegaFile);
    % omega_fixed.mat has variable 'omega_fixed': numYears x 6
    % columns: [year, All-N, High-N, High-U, High-E, All-E]
    omega_fixed_all = tmp.omega_fixed;
end

%% Build constraints (same for all years)
[A_ineq, B_ineq, Aeq, beq, lb, ub] = build_constraints(cfg);

%% Build initial guesses
candidates = build_initial_guesses(cfg, ub);
fprintf('Generated %d initial guess candidates\n', numel(candidates));

%% Set up parallel pool
pool = gcp('nocreate');
clust = parcluster('local');
maxW = clust.NumWorkers;
desiredW = min(cfg.numYears, maxW);
if isempty(pool)
    pool = parpool('local', desiredW);
end
pctRunOnAll addpath(pwd);
pctRunOnAll clear functions;
pctRunOnAll rng(0, 'twister');

%% Pre-allocate convergence tracking
ExitFlag   = nan(cfg.numYears, 1);
Fval       = nan(cfg.numYears, 1);
FirstOrder = nan(cfg.numYears, 1);
StepSz     = nan(cfg.numYears, 1);
Iterations = nan(cfg.numYears, 1);
Seconds    = nan(cfg.numYears, 1);

%% Main estimation loop
numYears = cfg.numYears;

parfor year = 1:numYears
    tStart = tic;
    fprintf('--- Year %d / %d ---\n', year, numYears);

    Data_full = Data_full_TimeSeries(:, year);

    % --- Build omega for this year ---
    if cfg.omegaEstimated
        omega_yr = [];  % will be extracted from x inside objective_fn
    else
        omeg1to5 = omega_fixed_all(year, 2:end)';
        % omeg1to5 layout: [All-N, High-N, High-U, High-E, All-E]
        omega_yr = zeros(cfg.ThetaTot, 1);
        if cfg.allE_estimated
            % 4 estimated types + 1 hardwired (All-N)
            omega_yr(1) = omeg1to5(2);  % High-N
            omega_yr(2) = 1 - omeg1to5(1) - omeg1to5(2) - omeg1to5(4) - omeg1to5(5);  % High-U as residual
            omega_yr(3) = omeg1to5(4);  % High-E
            omega_yr(4) = omeg1to5(5);  % All-E
            omega_yr(5) = omeg1to5(1);  % All-N (hardwired)
        else
            % 3 estimated types + 2 hardwired (All-E, All-N)
            omega_yr(1) = omeg1to5(2);  % High-N
            omega_yr(2) = omeg1to5(3);  % High-U
            omega_yr(3) = omeg1to5(4);  % High-E
            omega_yr(4) = omeg1to5(5);  % All-E (hardwired)
            omega_yr(5) = omeg1to5(1);  % All-N (hardwired)
        end
    end

    % --- Objective function handle ---
    f = @(x) objective_fn(x, cfg, Hstates, Data_full, omega_yr);

    % --- Pick best candidate ---
    bestx = candidates{1};
    bestf = inf;
    for c = 1:numel(candidates)
        try
            ftry = f(candidates{c});
        catch
            ftry = inf;
        end
        if ftry < bestf
            bestf = ftry;
            bestx = candidates{c};
        end
    end

    % --- Optimizer: interior-point (no ScaleProblem on first pass) ---
    options1 = optimoptions('fmincon', ...
        'Algorithm',             'interior-point', ...
        'Display',               'iter', ...
        'MaxIterations',         cfg.MaxIter, ...
        'MaxFunctionEvaluations',cfg.MaxEvals, ...
        'OptimalityTolerance',   cfg.OptTols, ...
        'StepTolerance',         1e-10, ...
        'ConstraintTolerance',   1e-9, ...
        'FiniteDifferenceType',  'forward', ...
        'FiniteDifferenceStepSize', 1e-5);

    [ests, SSR2, exitflag, outputF, ~, ~, hessian] = ...
        fmincon(f, bestx, A_ineq, B_ineq, Aeq, beq, lb, ub, [], options1);

    % --- SQP refinement if needed ---
    if exitflag == 2 || exitflag == 0 || ...
       (isfield(outputF,'firstorderopt') && outputF.firstorderopt > 1e-3)

        opts2 = optimoptions('fmincon', ...
            'Algorithm',             'sqp', ...
            'Display',               'iter', ...
            'MaxIterations',         cfg.MaxIter, ...
            'MaxFunctionEvaluations',cfg.MaxEvals, ...
            'OptimalityTolerance',   cfg.OptTols, ...
            'StepTolerance',         1e-10, ...
            'ConstraintTolerance',   1e-9, ...
            'FiniteDifferenceType',  'central', ...
            'FiniteDifferenceStepSize', 1e-4, ...
            'HessianApproximation',  'lbfgs', ...
            'ScaleProblem',          'obj-and-constr');

        [e2, f2, ef2, out2, ~, ~, hes2] = ...
            fmincon(f, ests, A_ineq, B_ineq, Aeq, beq, lb, ub, [], opts2);

        Nobs   = sum(Data_full);
        tolAbs = max(1, Nobs) * 1e-6;
        better = (f2 < SSR2 + tolAbs);
        kkt_ok = isfield(out2,'firstorderopt') && isfield(outputF,'firstorderopt') ...
                 && (out2.firstorderopt < outputF.firstorderopt);

        if better && (ef2 == 1 || kkt_ok)
            [ests, SSR2, exitflag, outputF, hessian] = deal(e2, f2, ef2, out2, hes2);
        end

        % If still not converged, try interior-point with central differences
        if exitflag == 2 || exitflag == 0
            opts3 = optimoptions('fmincon', ...
                'Algorithm',             'interior-point', ...
                'Display',               'iter', ...
                'MaxIterations',         cfg.MaxIter, ...
                'MaxFunctionEvaluations',cfg.MaxEvals, ...
                'OptimalityTolerance',   cfg.OptTols, ...
                'StepTolerance',         1e-10, ...
                'ConstraintTolerance',   1e-9, ...
                'FiniteDifferenceType',  'central', ...
                'FiniteDifferenceStepSize', 1e-5);

            [e3, f3, ef3, out3, ~, ~, hes3] = ...
                fmincon(f, ests, A_ineq, B_ineq, Aeq, beq, lb, ub, [], opts3);

            if f3 < SSR2
                [ests, SSR2, exitflag, outputF, hessian] = deal(e3, f3, ef3, out3, hes3);
            end
        end
    end

    fprintf('Year %d: exitflag=%d, fval=%.6f\n', year, exitflag, SSR2);

    % --- Extract results into named struct ---
    res = extract_results(ests, cfg, omega_yr, SSR2, exitflag);

    % --- Save per-year output ---
    outFile = fullfile(cfg.outputDir, sprintf('%s_%d.mat', cfg.outputPrefix, year));
    results_notReLabelled = res.results_notReLabelled;
    results = res.results_reLabelled;
    parsave(outFile, res, results_notReLabelled, results);

    % --- Track convergence ---
    ExitFlag(year)   = exitflag;
    Fval(year)       = SSR2;
    Iterations(year) = outputF.iterations;
    if isfield(outputF, 'firstorderopt')
        FirstOrder(year) = outputF.firstorderopt;
    end
    if isfield(outputF, 'stepsize')
        StepSz(year) = outputF.stepsize;
    end
    Seconds(year) = toc(tStart);
end

%% Save convergence summary
T = table(cfg.yearsList(:), ExitFlag, FirstOrder, StepSz, ...
    Iterations, Fval, Seconds, ...
    'VariableNames', {'Year','ExitFlag','FirstOrderOpt','StepSize', ...
                      'Iterations','Fval','Seconds'});
disp(T);
writetable(T, fullfile(cfg.outputDir, ...
    sprintf('convergence_%s_group%d.csv', cfg.mode, cfg.groupNum)));
save(fullfile(cfg.outputDir, ...
    sprintf('convergence_%s_group%d.mat', cfg.mode, cfg.groupNum)), 'T');

fprintf('\n=== Estimation complete ===\n');
fprintf('ExitFlag=1 in %d / %d years\n', sum(ExitFlag==1), cfg.numYears);

% =========================================================================
function parsave(outFile, res, results_notReLabelled, results) %#ok<INUSL>
% parsave  Save variables inside a parfor loop.
    save(outFile, 'res', 'results_notReLabelled', 'results');
end
