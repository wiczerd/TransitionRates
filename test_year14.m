% test_year14.m — Quick single-year test (year 14) for the rewrite
%
% Runs the estimation pipeline on year 14 only and compares against
% existing output files from the old code.

clc; clear all;
fprintf('=== Test: year 14, single-year estimation ===\n');

%% Load configuration (default: fixed_omega, allE_estimated = Step 3b)
config_Paper2;

%% Build Hstates (vectorized)
fprintf('Building Hstates...\n');
grids = cell(1, cfg.Nmonth);
[grids{:}] = ndgrid(1:cfg.Nstates);
Hstates = zeros(cfg.Nstatepath, cfg.Nmonth);
for k = 1:cfg.Nmonth
    Hstates(:,k) = grids{k}(:);
end
clear grids;

%% Load data
dataFile = fullfile(cfg.oldDataDir, cfg.inputFileNames{cfg.groupNum});
fprintf('Loading data from %s ...\n', dataFile);
Data_full_TimeSeries = xlsread(dataFile, cfg.dataRange);
fprintf('  Data size: %d x %d\n', size(Data_full_TimeSeries));

year = 14;
Data_full = Data_full_TimeSeries(:, year);
fprintf('Year index %d (calendar year %d)\n', year, cfg.yearsList(year));

%% Load fixed omegas
tmp = load(cfg.omegaFile);
omega_fixed_all = tmp.omega_fixed;
omeg1to5 = omega_fixed_all(year, 2:end)';

% Build omega vector
omega_yr = zeros(cfg.ThetaTot, 1);
if cfg.allE_estimated
    omega_yr(1) = omeg1to5(2);  % High-N
    omega_yr(2) = 1 - omeg1to5(1) - omeg1to5(2) - omeg1to5(4) - omeg1to5(5);  % High-U residual
    omega_yr(3) = omeg1to5(4);  % High-E
    omega_yr(4) = omeg1to5(5);  % All-E
    omega_yr(5) = omeg1to5(1);  % All-N
end
fprintf('Omegas: '); fprintf('%.4f ', omega_yr); fprintf('\n');

%% Build constraints
[A_ineq, B_ineq, Aeq, beq, lb, ub] = build_constraints(cfg);
fprintf('Constraints: %d ineq, %d eq, %d vars\n', size(A_ineq,1), size(Aeq,1), numel(lb));

%% Build initial guesses
candidates = build_initial_guesses(cfg, ub);
fprintf('Generated %d initial guess candidates\n', numel(candidates));

%% Evaluate candidates
f = @(x) objective_fn(x, cfg, Hstates, Data_full, omega_yr);

fprintf('\nEvaluating candidates:\n');
bestx = candidates{1}; bestf = inf;
for c = 1:numel(candidates)
    try
        ftry = f(candidates{c});
    catch ME
        fprintf('  Candidate %d: ERROR - %s\n', c, ME.message);
        ftry = inf;
    end
    fprintf('  Candidate %d: fval = %.6f\n', c, ftry);
    if ftry < bestf
        bestf = ftry; bestx = candidates{c};
    end
end
fprintf('Best candidate: fval = %.6f\n', bestf);

%% Run fmincon (interior-point)
fprintf('\n--- Running fmincon (interior-point) ---\n');
options1 = optimoptions('fmincon', ...
    'Algorithm',             'interior-point', ...
    'Display',               'iter', ...
    'MaxIterations',         cfg.MaxIter, ...
    'MaxFunctionEvaluations',cfg.MaxEvals, ...
    'OptimalityTolerance',   cfg.OptTols, ...
    'StepTolerance',         1e-10, ...
    'ConstraintTolerance',   1e-9, ...
    'FiniteDifferenceType',  'forward', ...
    'FiniteDifferenceStepSize', 1e-5, ...
    'ScaleProblem',          'obj-and-constr');

tic;
[ests, SSR2, exitflag, outputF, ~, ~, hessian] = ...
    fmincon(f, bestx, A_ineq, B_ineq, Aeq, beq, lb, ub, [], options1);
t1 = toc;
fprintf('Interior-point: exitflag=%d, fval=%.6f, time=%.1fs\n', exitflag, SSR2, t1);

%% SQP refinement if needed
if exitflag == 2 || exitflag == 0 || ...
   (isfield(outputF,'firstorderopt') && outputF.firstorderopt > 1e-3)
    fprintf('\n--- Running SQP refinement ---\n');
    opts2 = optimoptions('fmincon', ...
        'Algorithm','sqp', ...
        'Display','iter', ...
        'MaxIterations',         cfg.MaxIter, ...
        'MaxFunctionEvaluations',cfg.MaxEvals, ...
        'OptimalityTolerance',   cfg.OptTols, ...
        'StepTolerance',         1e-10, ...
        'ConstraintTolerance',   1e-9, ...
        'FiniteDifferenceType','central', ...
        'FiniteDifferenceStepSize',1e-4, ...
        'HessianApproximation','lbfgs', ...
        'ScaleProblem','obj-and-constr');

    tic;
    [e2, f2, ef2, out2, ~, ~, hes2] = ...
        fmincon(f, ests, A_ineq, B_ineq, Aeq, beq, lb, ub, [], opts2);
    t2 = toc;
    fprintf('SQP: exitflag=%d, fval=%.6f, time=%.1fs\n', ef2, f2, t2);

    Nobs   = sum(Data_full);
    tolAbs = max(1, Nobs) * 1e-6;
    better = (f2 < SSR2 + tolAbs);
    kkt_ok = isfield(out2,'firstorderopt') && isfield(outputF,'firstorderopt') ...
             && (out2.firstorderopt < outputF.firstorderopt);
    if better && (ef2 == 1 || kkt_ok)
        [ests, SSR2, exitflag, outputF, hessian] = deal(e2, f2, ef2, out2, hes2);
        fprintf('Accepted SQP solution\n');
    end
end

%% Extract results
res = extract_results(ests, cfg, omega_yr, SSR2, exitflag);

fprintf('\n=== Results for year %d ===\n', year);
fprintf('Exit flag: %d\n', res.exitflag);
fprintf('Objective: %.6f\n', res.SSR2);
fprintf('Omega: '); fprintf('%.4f ', res.omega); fprintf('\n');
tNames = {'High-N','High-U','High-E','All-E'};
for th = 1:cfg.Theta
    fprintf('\nType %d (%s):\n', th, tNames{th});
    fprintf('  Ergodic: [%.4f %.4f %.4f %.4f]\n', res.ergVec{th});
    fprintf('  U-rate:  %.4f\n', res.urate(th));
    fprintf('  TM:\n');
    disp(res.tm{th});
end
fprintf('Aggregate ergodic: [%.4f %.4f %.4f %.4f]\n', res.ergAgg);

%% Save test output
outFile = fullfile(cfg.outputDir, sprintf('%s_%d.mat', cfg.outputPrefix, year));
results_notReLabelled = res.results_notReLabelled;
results = res.results_reLabelled;
save(outFile, 'res', 'results_notReLabelled', 'results');
fprintf('\nSaved to %s\n', outFile);

%% Compare with old output if available
oldFile = fullfile(cfg.oldOutputDir, sprintf('outputFO_notReLabelled_AllE_V3_%d.mat', year));
if isfile(oldFile)
    fprintf('\n=== Comparison with old code ===\n');
    old = load(oldFile);
    old_vec = old.results_notReLabelled;
    new_vec = res.results_notReLabelled;
    fprintf('Old output length: %d\n', numel(old_vec));
    fprintf('New output length: %d\n', numel(new_vec));
    if numel(old_vec) == numel(new_vec)
        maxdiff = max(abs(old_vec - new_vec));
        fprintf('Max absolute difference: %.6e\n', maxdiff);
        if maxdiff < 1e-3
            fprintf('MATCH: Results agree within 1e-3\n');
        else
            fprintf('DIVERGENT: Results differ (expected if different starting point)\n');
            fprintf('Old fval (SSR2): %.6f\n', old_vec(13));
            fprintf('New fval (SSR2): %.6f\n', new_vec(13));
        end
    else
        fprintf('Different vector lengths - comparing fval only\n');
        % In old format, SSR2 is at position 13 of the packed column
        % (12 pVec entries for type 1 + 1 SSR2 entry)
        fprintf('Old fval: check manually\n');
        fprintf('New fval: %.6f\n', SSR2);
    end
else
    fprintf('\nNo old output file found for comparison at:\n  %s\n', oldFile);
end

fprintf('\n=== Test complete ===\n');
