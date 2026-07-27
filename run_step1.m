% run_step1.m — Step 1: Free omegas estimation, all years
%
% Estimates per year: type shares (omega) + transition matrices
% 5 types = 3 mover (High-N, High-U, High-E) with full 4x4 matrices
%         + 2 polar (All-N, All-E) hardwired
%
% Toggle: config_Paper2.m → cfg.mode = 'free_omega'
%
% Output: Step1_FreeOmega/ subfolder with per-year .mat, Excel, figures

clc; clear all;
tTotal = tic;

%% Load configuration — force free_omega mode
config_Paper2;
cfg.mode = 'free_omega';
cfg.Theta    = 3;
cfg.ThetaTot = 5;
cfg.omegaEstimated = true;
cfg.allE_estimated = false;
cfg.useRegularization = false;  % no regularization for Step 1
cfg.outputPrefix = 'step1';
cfg.outputDir = fullfile(cfg.dataDir, 'Step1_FreeOmega');

fprintf('=== Step 1: Free Omegas Estimation ===\n');
fprintf('Theta=%d mover types + 2 polar = %d total, %d years\n', ...
    cfg.Theta, cfg.ThetaTot, cfg.numYears);

%% Build Hstates (vectorized with ndgrid)
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
fprintf('  Data: %d paths x %d years\n', size(Data_full_TimeSeries));

%% Build constraints
[A_ineq, B_ineq, Aeq, beq, lb, ub] = build_constraints(cfg);
fprintf('Parameters: %d (36 transition + 5 omegas)\n', numel(lb));

%% Build initial guesses
candidates = build_initial_guesses(cfg, ub);
fprintf('Initial guess candidates: %d\n', numel(candidates));

%% Set up parallel pool
pool = gcp('nocreate');
% Query the local cluster's max workers (respects MATLAB license/settings)
clust = parcluster('local');
maxW = clust.NumWorkers;
desiredWorkers = min(cfg.numYears, maxW);
if isempty(pool)
    pool = parpool('local', desiredWorkers);
elseif pool.NumWorkers ~= desiredWorkers
    delete(pool);
    pool = parpool('local', desiredWorkers);
end
fprintf('Parallel pool: %d workers\n', pool.NumWorkers);
pctRunOnAll addpath(fileparts(mfilename('fullpath')));
pctRunOnAll clear functions;
pctRunOnAll rng(0, 'twister');

%% Pre-allocate
numYears = cfg.numYears;
ExitFlag   = nan(numYears, 1);
Fval       = nan(numYears, 1);
FirstOrder = nan(numYears, 1);
Iterations = nan(numYears, 1);
Seconds    = nan(numYears, 1);

%% Main estimation loop
fprintf('\n--- Starting estimation for %d years ---\n', numYears);

parfor year = 1:numYears
    tStart = tic;
    Data_full = Data_full_TimeSeries(:, year);
    omega_yr = [];  % extracted from x in objective_fn

    f = @(x) objective_fn(x, cfg, Hstates, Data_full, omega_yr);

    % Pick best candidate
    bestx = candidates{1}; bestf = inf;
    for c = 1:numel(candidates)
        try, ftry = f(candidates{c}); catch, ftry = inf; end
        if ftry < bestf, bestf = ftry; bestx = candidates{c}; end
    end

    % Interior-point (no ScaleProblem — matches old code for stable omega gradient)
    opts1 = optimoptions('fmincon', ...
        'Algorithm','interior-point', 'Display','off', ...
        'MaxIterations',cfg.MaxIter, 'MaxFunctionEvaluations',cfg.MaxEvals, ...
        'OptimalityTolerance',cfg.OptTols, 'StepTolerance',1e-10, ...
        'ConstraintTolerance',1e-9, 'FiniteDifferenceType','forward', ...
        'FiniteDifferenceStepSize',1e-5);

    [ests, SSR2, exitflag, outputF, ~, ~, hessian] = ...
        fmincon(f, bestx, A_ineq, B_ineq, Aeq, beq, lb, ub, [], opts1);

    % SQP refinement
    if exitflag ~= 1 || (isfield(outputF,'firstorderopt') && outputF.firstorderopt > 1e-3)
        opts2 = optimoptions('fmincon', ...
            'Algorithm','sqp', 'Display','off', ...
            'MaxIterations',cfg.MaxIter, 'MaxFunctionEvaluations',cfg.MaxEvals, ...
            'OptimalityTolerance',cfg.OptTols, 'StepTolerance',1e-10, ...
            'ConstraintTolerance',1e-9, 'FiniteDifferenceType','central', ...
            'FiniteDifferenceStepSize',1e-4, 'HessianApproximation','lbfgs', ...
            'ScaleProblem','obj-and-constr');
        [e2,f2,ef2,out2,~,~,hes2] = fmincon(f,ests,A_ineq,B_ineq,Aeq,beq,lb,ub,[],opts2);
        if f2 <= SSR2 && (ef2 == 1 || f2 < SSR2 - 1e-6)
            [ests,SSR2,exitflag,outputF,hessian] = deal(e2,f2,ef2,out2,hes2);
        end
    end

    % Interior-point central as last resort
    if exitflag ~= 1
        opts3 = optimoptions('fmincon', ...
            'Algorithm','interior-point', 'Display','off', ...
            'MaxIterations',cfg.MaxIter, 'MaxFunctionEvaluations',cfg.MaxEvals, ...
            'OptimalityTolerance',cfg.OptTols, 'StepTolerance',1e-10, ...
            'ConstraintTolerance',1e-9, 'FiniteDifferenceType','central', ...
            'FiniteDifferenceStepSize',1e-5);
        [e3,f3,ef3,out3,~,~,hes3] = fmincon(f,ests,A_ineq,B_ineq,Aeq,beq,lb,ub,[],opts3);
        if f3 < SSR2
            [ests,SSR2,exitflag,outputF,hessian] = deal(e3,f3,ef3,out3,hes3);
        end
    end

    % Extract results
    res = extract_results(ests, cfg, omega_yr, SSR2, exitflag);

    % Save
    outFile = fullfile(cfg.outputDir, sprintf('%s_%d.mat', cfg.outputPrefix, year));
    parsave_step1(outFile, res);

    % Track convergence
    ExitFlag(year) = exitflag;
    Fval(year) = SSR2;
    if isfield(outputF,'iterations'), Iterations(year) = outputF.iterations; end
    if isfield(outputF,'firstorderopt'), FirstOrder(year) = outputF.firstorderopt; end
    Seconds(year) = toc(tStart);

    fprintf('Year %d (%d): ef=%d, fval=%.4f, %.0fs\n', ...
        year, cfg.yearsList(year), exitflag, SSR2, Seconds(year));
end

%% Convergence summary
T = table(cfg.yearsList(:), ExitFlag, FirstOrder, Iterations, Fval, Seconds, ...
    'VariableNames', {'Year','ExitFlag','FirstOrderOpt','Iterations','Fval','Seconds'});
disp(T);
writetable(T, fullfile(cfg.outputDir, 'convergence_step1.csv'));
save(fullfile(cfg.outputDir, 'convergence_step1.mat'), 'T', 'cfg');

fprintf('\n=== Step 1 complete ===\n');
fprintf('ExitFlag=1: %d / %d years\n', sum(ExitFlag==1), numYears);
fprintf('Total time: %.1f minutes\n', toc(tTotal)/60);

% After estimation, run export and figures
fprintf('\nExporting to Excel...\n');
export_step1;
fprintf('Generating figures...\n');
figures_step1;

function parsave_step1(outFile, res)
    save(outFile, 'res');
end
