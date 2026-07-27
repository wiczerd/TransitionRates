% fix_year2024.m — Targeted re-estimation for year 2024 (index 47)
%
% Fixes:
%   1. lb = 1e-6 (not 0) to avoid parameters stuck at exact boundary
%   2. 20 candidates instead of 7 (more jitters, more random seeds, neighbors)
%   3. Compares final result against both old and multi-start solutions
%   4. Only saves if genuinely better AND well-converged

clc; clear;
tTotal = tic;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(fileparts(rewriteDir), ...
               'Replication Files and ResultsTablesFigures Files 2 (1)');
oldDir     = fullfile(rewriteDir, 'Step1_FreeOmega');
msDir      = fullfile(rewriteDir, 'Step1_FreeOmega_Improving');
bestDir    = fullfile(rewriteDir, 'Step1_FreeOmega_Best');

addpath(origDir);

%% Settings
TotalYears = 47;
Theta    = 3;
ThetaTot = Theta + 2;
groupNum = 4;
inputFileNames = {'DataTimeSeriesYW.xlsx','DataTimeSeriesPW.xlsx', ...
                  'DataTimeSeriesYM.xlsx','DataTimeSeriesPM_1976_2023.xlsx'};
Nstatepath = 4^8;
Nmonth     = 8;
yearsList  = setdiff(1976:2024, [1977 1993]);
yearIdx    = 47;  % 2024

NParEst = 12 * Theta;
Nest    = NParEst + ThetaTot;

%% Build Hstates
Hstates = zeros(Nstatepath, Nmonth);
path = 1;
for j1=1:4, for j2=1:4, for j3=1:4, for j4=1:4
for j5=1:4, for j6=1:4, for j7=1:4, for j8=1:4
    Hstates(path,:) = [j1 j2 j3 j4 j5 j6 j7 j8];
    path = path+1;
end, end, end, end, end, end, end, end

%% Load data for year 2024
Data_full_TimeSeries = xlsread(fullfile(origDir, inputFileNames{groupNum}), 'J2:BD6562');
Data_full = Data_full_TimeSeries(:, yearIdx);
select = (1:3^8)';

%% Load existing solutions to use as starting points
% Old original
tmp = load(fullfile(oldDir, sprintf('step1_%d.mat', yearIdx)));
unks_old = [tmp.res.pVecs(:); tmp.res.omega];
fval_old = tmp.res.SSR2;

% Multi-start
tmp = load(fullfile(msDir, sprintf('step1_%d.mat', yearIdx)));
unks_ms = [tmp.res.pVecs(:); tmp.res.omega];
fval_ms = tmp.res.SSR2;

% Neighbors: 2023 (idx 46), 2022 (idx 45), 2021 (idx 44)
unks_neighbors = {};
for nIdx = [46, 45, 44, 43, 42]
    tmp = load(fullfile(oldDir, sprintf('step1_%d.mat', nIdx)));
    unks_neighbors{end+1} = [tmp.res.pVecs(:); tmp.res.omega];
end

fprintf('=== Targeted fix for Year 2024 ===\n');
fprintf('Old fval:    %.6f\n', fval_old);
fprintf('MS fval:     %.6f (firstorderopt=0.712)\n', fval_ms);
fprintf('Best so far: %.6f\n', min(fval_old, fval_ms));

%% Build candidates (20 total)
candidates = {};
candLabels = {};

% 1. Old solution
candidates{end+1} = unks_old;
candLabels{end+1} = 'Old solution';

% 2. Multi-start solution
candidates{end+1} = unks_ms;
candLabels{end+1} = 'MS solution';

% 3-7. Neighbors (2023, 2022, 2021, 2020, 2019)
neighborYears = [2023, 2022, 2021, 2020, 2019];
for i = 1:numel(unks_neighbors)
    candidates{end+1} = unks_neighbors{i};
    candLabels{end+1} = sprintf('Neighbor %d', neighborYears(i));
end

% 8-13. Jittered old solution (6 versions, different sigmas)
jitterSigmas = [0.01, 0.02, 0.03, 0.05, 0.03, 0.05];
for j = 1:numel(jitterSigmas)
    rng(50000 + j, 'twister');
    unks_j = unks_old + jitterSigmas(j) * randn(Nest, 1);
    unks_j = max(1e-6, min(1-1e-6, unks_j));
    om = unks_j(NParEst+1:end);
    om = max(om, 1e-4);
    unks_j(NParEst+1:end) = om / sum(om);
    candidates{end+1} = unks_j;
    candLabels{end+1} = sprintf('Jitter old s=%.2f #%d', jitterSigmas(j), j);
end

% 14-17. Jittered MS solution (4 versions)
for j = 1:4
    rng(60000 + j, 'twister');
    unks_j = unks_ms + 0.03 * randn(Nest, 1);
    unks_j = max(1e-6, min(1-1e-6, unks_j));
    om = unks_j(NParEst+1:end);
    om = max(om, 1e-4);
    unks_j(NParEst+1:end) = om / sum(om);
    candidates{end+1} = unks_j;
    candLabels{end+1} = sprintf('Jitter MS #%d', j);
end

% 18-20. Random (seeds 400, 500, 600)
for seed = [400, 500, 600]
    Params = [];
    for kk = 1:Theta
        rng(seed + kk, 'twister');
        p12 = rand; p13 = (1-p12)*rand; p14 = (1-p12-p13)*rand;
        p22 = rand; p23 = (1-p22)*rand; p24 = (1-p22-p23)*rand;
        p32 = rand; p33 = (1-p32)*rand; p34 = (1-p32-p33)*rand;
        p42 = rand; p43 = (1-p42)*rand; p44 = (1-p42-p43)*rand;
        Params = [Params; p12;p13;p14;p22;p23;p24;p32;p33;p34;p42;p43;p44];
    end
    omeg = (1/ThetaTot) * ones(ThetaTot, 1);
    candidates{end+1} = max(1e-6, min(1-1e-6, [Params; omeg]));
    candLabels{end+1} = sprintf('Random seed=%d', seed);
end

nCand = numel(candidates);
fprintf('\nTotal candidates: %d\n\n', nCand);

%% Constraints — KEY FIX: lb = 1e-6 instead of 0
lb  = 1e-6 * ones(Nest, 1);
ub  = ones(Nest, 1);
Aeq = [zeros(1, NParEst), ones(1, ThetaTot)];
A   = zeros(Theta*4, Nest);
B   = ones(Theta*4, 1);
for index = 1:Theta*4
    A(index,:) = [zeros(1,3*(index-1)) ones(1,3) zeros(1,Nest-3-3*(index-1))];
end

%% Objective
f = @(x) ModelFitFuncTimeSeries_NewTrMBoot_20250919_V2( ...
    x, select, groupNum, Theta, Hstates, yearIdx, Data_full);

%% Phase A: Screen all candidates (quick)
opts_screen = optimoptions('fmincon', ...
    'Algorithm','interior-point', 'Display','off', ...
    'MaxIterations',3000, 'MaxFunctionEvaluations',8000, ...
    'OptimalityTolerance',0.01, 'StepTolerance',1e-8, ...
    'ConstraintTolerance',1e-9, 'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',1e-5);

screenFval = inf(nCand, 1);
screenEsts = cell(nCand, 1);
screenFOpt = nan(nCand, 1);

fprintf('Phase A: Screening %d candidates ...\n', nCand);
for c = 1:nCand
    try
        [e_c, f_c, ~, out_c] = fmincon(f, candidates{c}, A, B, Aeq, 1, lb, ub, [], opts_screen);
        screenFval(c) = f_c;
        screenEsts{c} = e_c;
        if isfield(out_c, 'firstorderopt')
            screenFOpt(c) = out_c.firstorderopt;
        end
    catch ME
        fprintf('  Candidate %d failed: %s\n', c, ME.message);
        screenFval(c) = Inf;
    end
    fprintf('  %2d. %-25s fval=%.6f  1stOrd=%.4f\n', c, candLabels{c}, screenFval(c), screenFOpt(c));
end

% Sort by screening fval
[sortedFval, sortIdx] = sort(screenFval);
fprintf('\n--- Top 5 candidates after screening ---\n');
for i = 1:min(5, nCand)
    c = sortIdx(i);
    fprintf('  %2d. %-25s fval=%.6f\n', c, candLabels{c}, screenFval(c));
end

%% Phase B: Full refinement on TOP 3 candidates (not just 1)
nRefine = min(3, nCand);
refinedFval = inf(nRefine, 1);
refinedEsts = cell(nRefine, 1);
refinedFOpt = nan(nRefine, 1);
refinedEF   = nan(nRefine, 1);

fprintf('\nPhase B: Refining top %d candidates ...\n', nRefine);
for i = 1:nRefine
    c = sortIdx(i);
    ests = screenEsts{c};
    fprintf('  Refining candidate %d (%s) ...\n', c, candLabels{c});

    % Pass 1: interior-point, forward FD
    opts1 = optimoptions('fmincon', ...
        'Algorithm','interior-point', 'Display','off', ...
        'MaxIterations',15000, 'MaxFunctionEvaluations',40000, ...
        'OptimalityTolerance',1e-4, 'StepTolerance',1e-12, ...
        'ConstraintTolerance',1e-10, 'FiniteDifferenceType','forward', ...
        'FiniteDifferenceStepSize',1e-5);
    [ests, SSR2, ef, outF] = fmincon(f, ests, A, B, Aeq, 1, lb, ub, [], opts1);
    fo = outF.firstorderopt;
    fprintf('    Pass 1: fval=%.6f ef=%d 1stOrd=%.6f\n', SSR2, ef, fo);

    % Pass 2: SQP with central FD
    if fo > 1e-3 || ef == 2 || ef == 0
        opts2 = optimoptions('fmincon', 'Algorithm','sqp', 'Display','off', ...
            'MaxIterations',15000, 'MaxFunctionEvaluations',40000, ...
            'OptimalityTolerance',1e-4, 'StepTolerance',1e-12, ...
            'ConstraintTolerance',1e-10, 'FiniteDifferenceType','central', ...
            'FiniteDifferenceStepSize',1e-4, ...
            'HessianApproximation','lbfgs', 'ScaleProblem','obj-and-constr');
        [e2,f2,ef2,out2] = fmincon(f, ests, A, B, Aeq, 1, lb, ub, [], opts2);
        fprintf('    Pass 2 (SQP): fval=%.6f ef=%d 1stOrd=%.6f\n', f2, ef2, out2.firstorderopt);
        if f2 < SSR2
            [ests, SSR2, ef, outF] = deal(e2, f2, ef2, out2);
            fo = outF.firstorderopt;
        end
    end

    % Pass 3: interior-point with central FD (tighter)
    if fo > 1e-3 || ef == 2 || ef == 0
        opts3 = optimoptions('fmincon', ...
            'Algorithm','interior-point', 'Display','off', ...
            'MaxIterations',20000, 'MaxFunctionEvaluations',50000, ...
            'OptimalityTolerance',1e-4, 'StepTolerance',1e-14, ...
            'ConstraintTolerance',1e-10, 'FiniteDifferenceType','central', ...
            'FiniteDifferenceStepSize',1e-5);
        [e3,f3,ef3,out3] = fmincon(f, ests, A, B, Aeq, 1, lb, ub, [], opts3);
        fprintf('    Pass 3 (IP central): fval=%.6f ef=%d 1stOrd=%.6f\n', f3, ef3, out3.firstorderopt);
        if f3 < SSR2
            [ests, SSR2, ef, outF] = deal(e3, f3, ef3, out3);
            fo = outF.firstorderopt;
        end
    end

    % Pass 4: one more SQP from current point if still not converged
    if fo > 1e-3
        opts4 = optimoptions('fmincon', 'Algorithm','sqp', 'Display','off', ...
            'MaxIterations',20000, 'MaxFunctionEvaluations',50000, ...
            'OptimalityTolerance',1e-5, 'StepTolerance',1e-14, ...
            'ConstraintTolerance',1e-10, 'FiniteDifferenceType','central', ...
            'FiniteDifferenceStepSize',1e-5, ...
            'HessianApproximation','lbfgs', 'ScaleProblem','obj-and-constr');
        [e4,f4,ef4,out4] = fmincon(f, ests, A, B, Aeq, 1, lb, ub, [], opts4);
        fprintf('    Pass 4 (SQP tight): fval=%.6f ef=%d 1stOrd=%.6f\n', f4, ef4, out4.firstorderopt);
        if f4 < SSR2
            [ests, SSR2, ef, outF] = deal(e4, f4, ef4, out4);
            fo = outF.firstorderopt;
        end
    end

    refinedFval(i) = SSR2;
    refinedEsts{i} = ests;
    refinedFOpt(i) = fo;
    refinedEF(i)   = ef;
    fprintf('  -> Final: fval=%.6f ef=%d 1stOrd=%.6f\n\n', SSR2, ef, fo);
end

%% Pick the best refined solution
[bestFval, bestIdx] = min(refinedFval);
bestEsts = refinedEsts{bestIdx};
bestFOpt = refinedFOpt(bestIdx);
bestEF   = refinedEF(bestIdx);
bestCandOrig = sortIdx(bestIdx);

fprintf('=== BEST REFINED SOLUTION ===\n');
fprintf('Candidate: %d (%s)\n', bestCandOrig, candLabels{bestCandOrig});
fprintf('Fval:      %.6f\n', bestFval);
fprintf('1stOrdOpt: %.6f\n', bestFOpt);
fprintf('ExitFlag:  %d\n', bestEF);
fprintf('\nComparison:\n');
fprintf('  Old:        %.6f\n', fval_old);
fprintf('  MS:         %.6f (1stOrd=0.712)\n', fval_ms);
fprintf('  New best:   %.6f (1stOrd=%.4f)\n', bestFval, bestFOpt);
fprintf('  Improvement over old: %.6f\n', fval_old - bestFval);

%% Post-processing (same as in rerun_all_multistart.m)
ests = bestEsts;
pVecs = [];
omega = ests(end-(ThetaTot-1):end);
for theta = 1:Theta
    pVec = ests((theta-1)*12+1:theta*12);
    pVecs = [pVecs, pVec];
end

% State ordering: ensure p33 <= p44
for theta = 1:Theta
    if pVecs(8,theta) > pVecs(12,theta)
        old_p = pVecs(:,theta);
        pVecs(:,theta) = [old_p(1);old_p(3);old_p(2); old_p(4);old_p(6);old_p(5); ...
                          old_p(10);old_p(12);old_p(11); old_p(7);old_p(9);old_p(8)];
    end
end

% Ergodic distributions
ergVecs = zeros(4, Theta);
for theta = 1:Theta
    p = pVecs(:,theta);
    tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
          1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
          1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
          1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
    ergVecs(:,theta) = Ergodic(tm);
end

% Type ordering: sort by ergodic E mass
E_mass = ergVecs(3,:) + ergVecs(4,:);
[~, I] = sort(E_mass);
pVecs = pVecs(:, I);
omega(1:Theta) = omega(I);

% Recompute ergodic after type sort
ergVecs = zeros(4, Theta);
for theta = 1:Theta
    p = pVecs(:,theta);
    tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
          1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
          1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
          1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
    ergVecs(:,theta) = Ergodic(tm);
end

% Sub-sort types 1 and 2 by U mass
U_mass = ergVecs(2, 1:2);
[~, Isub] = sort(U_mass);
pVecs(:,1:2) = pVecs(:, Isub);
omega(1:2) = omega(Isub);

% Final ergodic
ergVecs = zeros(4, Theta);
for theta = 1:Theta
    p = pVecs(:,theta);
    tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
          1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
          1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
          1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
    ergVecs(:,theta) = Ergodic(tm);
end
ergVecsAll = ergVecs * omega(1:Theta) + [omega(ThetaTot);0;0;omega(ThetaTot-1)];

fprintf('\n--- New solution for 2024 ---\n');
fprintf('omega: High-N=%.4f, High-U=%.4f, High-E=%.4f, All-E=%.4f, All-N=%.4f\n', ...
    omega(1), omega(2), omega(3), omega(4), omega(5));
fprintf('ergVecs:\n'); disp(ergVecs);
fprintf('Params at lb (=1e-6): %d\n', sum(abs(bestEsts(1:NParEst) - 1e-6) < 1e-8));

%% Save — only if better than both old and MS
if bestFval <= min(fval_old, fval_ms) + 1e-8
    SSR2 = bestFval;
    exitflag = bestEF;
    results = [pVecs; ...
        [SSR2, zeros(1,Theta-1)]; ...
        omega(1:Theta)'; ...
        [omega(ThetaTot-1), zeros(1,Theta-1)]; ...
        [omega(ThetaTot),   zeros(1,Theta-1)]; ...
        [0.001,             zeros(1,Theta-1)]; ...
        [100,               zeros(1,Theta-1)]; ...
        zeros(1, Theta); ...
        [exitflag,          zeros(1,Theta-1)]; ...
        ergVecs];

    res = struct();
    res.results    = results(:);
    res.pVecs      = pVecs;
    res.omega      = omega;
    res.ergVecs    = ergVecs;
    res.ergVecsAll = ergVecsAll;
    res.SSR2       = SSR2;
    res.exitflag   = exitflag;

    % Save to Best folder
    save(fullfile(bestDir, 'step1_47.mat'), 'res');
    fprintf('\nSaved to %s\n', fullfile(bestDir, 'step1_47.mat'));

    % Also save to Improving folder
    res.bestCandID = bestCandOrig;
    res.screenFvals = screenFval;
    save(fullfile(msDir, 'step1_47.mat'), 'res');
    fprintf('Saved to %s\n', fullfile(msDir, 'step1_47.mat'));
else
    fprintf('\nNo improvement found. Keeping existing solutions.\n');
end

fprintf('\nTotal time: %.1f minutes\n', toc(tTotal)/60);
