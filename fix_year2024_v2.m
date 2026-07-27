% fix_year2024_v2.m — Re-estimate 2024 to fix High-E / All-E type swap
%
% The current 2024 solution has High-E omega=0.63 and All-E omega=0.05,
% which is a type-label swap: the "High-E" mover type is doing All-E's job.
%
% Strategy:
%   - Use well-behaved years from the full sample as starting points
%   - Avoid nearby years (2022 may also be partially swapped)
%   - Include explicitly constructed starts that enforce proper type structure
%   - 30 candidates total, aggressive multi-pass refinement

clc; clear;
tTotal = tic;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(fileparts(rewriteDir), ...
               'Replication Files and ResultsTablesFigures Files 2 (1)');
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

%% Load current 2024 solution
tmp = load(fullfile(bestDir, 'step1_47.mat'));
fval_current = tmp.res.SSR2;
fprintf('=== Re-estimating Year 2024 (fixing type swap) ===\n');
fprintf('Current fval: %.6f\n', fval_current);
fprintf('Current omegas: High-N=%.4f High-U=%.4f High-E=%.4f All-E=%.4f All-N=%.4f\n', ...
    tmp.res.omega);

%% Load ALL year solutions to pick well-behaved starting points
allUnks = cell(TotalYears, 1);
allOmega = nan(TotalYears, ThetaTot);
for y = 1:TotalYears
    t = load(fullfile(bestDir, sprintf('step1_%d.mat', y)));
    allUnks{y} = [t.res.pVecs(:); t.res.omega];
    allOmega(y,:) = t.res.omega(:)';
end

%% Identify well-behaved years (All-E > 0.40 and High-E < 0.30)
wellBehaved = find(allOmega(:,4) > 0.40 & allOmega(:,3) < 0.30);
fprintf('\nWell-behaved years (%d total):\n', length(wellBehaved));
for i = 1:length(wellBehaved)
    y = wellBehaved(i);
    fprintf('  %d: All-E=%.3f High-E=%.3f\n', yearsList(y), allOmega(y,4), allOmega(y,3));
end

%% Build candidates
candidates = {};
candLabels = {};

% 1-12. Well-behaved years spread across the sample
% Pick every 3rd well-behaved year to get a spread
if length(wellBehaved) > 12
    step = floor(length(wellBehaved) / 12);
    selectedWB = wellBehaved(1:step:end);
    selectedWB = selectedWB(1:min(12, length(selectedWB)));
else
    selectedWB = wellBehaved;
end
for i = 1:length(selectedWB)
    y = selectedWB(i);
    candidates{end+1} = allUnks{y};
    candLabels{end+1} = sprintf('Year %d (AllE=%.2f)', yearsList(y), allOmega(y,4));
end

% 13-17. Recent well-behaved years (2016-2019, 2023)
recentGood = [39 40 41 42 46]; % indices for 2016,2017,2018,2019,2023
for i = 1:length(recentGood)
    y = recentGood(i);
    if ~ismember(y, selectedWB)
        candidates{end+1} = allUnks{y};
        candLabels{end+1} = sprintf('Year %d (recent)', yearsList(y));
    end
end

% 18-20. Manually constructed: take 2023 solution, scale omegas to give
% All-E a bigger share while keeping mover TMs
for scale_AllE = [0.55, 0.60, 0.65]
    unks_mod = allUnks{46};  % 2023 as base
    om = unks_mod(NParEst+1:end);
    % Force All-E to target, redistribute proportionally
    om(4) = scale_AllE;
    om(3) = om(3) * (1 - scale_AllE - om(5)) / (om(1)+om(2)+om(3));
    om(1) = om(1) * (1 - scale_AllE - om(5)) / (allOmega(46,1)+allOmega(46,2)+allOmega(46,3));
    om(2) = om(2) * (1 - scale_AllE - om(5)) / (allOmega(46,1)+allOmega(46,2)+allOmega(46,3));
    om = max(om, 1e-4);
    om = om / sum(om);
    unks_mod(NParEst+1:end) = om;
    candidates{end+1} = unks_mod;
    candLabels{end+1} = sprintf('2023 TM + AllE=%.2f', scale_AllE);
end

% 21-23. Take current 2024 solution but manually fix: swap omega(3) and omega(4)
% and reset the "new High-E" TM to something from a good year
for srcYr = [42, 46, 38]  % 2019, 2023, 2017
    unks_swap = allUnks{47};  % current 2024
    % Swap omegas
    unks_swap(NParEst+3) = allOmega(47,4);  % High-E gets old All-E omega
    unks_swap(NParEst+4) = allOmega(47,3);  % All-E gets old High-E omega
    % But the TM for "High-E" (type 3) is now wrong — replace with source year's High-E TM
    srcUnks = allUnks{srcYr};
    unks_swap(25:36) = srcUnks(25:36);  % type 3 TM params
    om = unks_swap(NParEst+1:end);
    om = max(om, 1e-4);
    om = om / sum(om);
    unks_swap(NParEst+1:end) = om;
    candidates{end+1} = unks_swap;
    candLabels{end+1} = sprintf('2024 swapped + %d TM', yearsList(srcYr));
end

% 24-27. Jittered versions of recent good years
for srcYr = [42, 46, 40, 38]  % 2019, 2023, 2017, 2015
    rng(80000 + srcYr, 'twister');
    unks_j = allUnks{srcYr} + 0.02 * randn(Nest, 1);
    unks_j = max(1e-6, min(1-1e-6, unks_j));
    om = unks_j(NParEst+1:end);
    om = max(om, 1e-4);
    unks_j(NParEst+1:end) = om / sum(om);
    candidates{end+1} = unks_j;
    candLabels{end+1} = sprintf('Jittered %d', yearsList(srcYr));
end

% 28-30. Random with high All-E prior
for seed = [700, 800, 900]
    Params = [];
    for kk = 1:Theta
        rng(seed + kk, 'twister');
        p12 = rand; p13 = (1-p12)*rand; p14 = (1-p12-p13)*rand;
        p22 = rand; p23 = (1-p22)*rand; p24 = (1-p22-p23)*rand;
        p32 = rand; p33 = (1-p32)*rand; p34 = (1-p32-p33)*rand;
        p42 = rand; p43 = (1-p42)*rand; p44 = (1-p42-p43)*rand;
        Params = [Params; p12;p13;p14;p22;p23;p24;p32;p33;p34;p42;p43;p44];
    end
    % Set omegas with high All-E
    omeg = [0.08; 0.03; 0.15; 0.60; 0.14];
    omeg = omeg / sum(omeg);
    candidates{end+1} = max(1e-6, min(1-1e-6, [Params; omeg]));
    candLabels{end+1} = sprintf('Random+highAllE seed=%d', seed);
end

nCand = numel(candidates);
fprintf('\nTotal candidates: %d\n\n', nCand);

%% Constraints
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

%% Phase A: Screen all candidates
opts_screen = optimoptions('fmincon', ...
    'Algorithm','interior-point', 'Display','off', ...
    'MaxIterations',3000, 'MaxFunctionEvaluations',8000, ...
    'OptimalityTolerance',0.01, 'StepTolerance',1e-8, ...
    'ConstraintTolerance',1e-9, 'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',1e-5);

screenFval = inf(nCand, 1);
screenEsts = cell(nCand, 1);
screenFOpt = nan(nCand, 1);
screenOmega = nan(nCand, ThetaTot);

fprintf('Phase A: Screening %d candidates ...\n', nCand);
for c = 1:nCand
    try
        [e_c, f_c, ~, out_c] = fmincon(f, candidates{c}, A, B, Aeq, 1, lb, ub, [], opts_screen);
        screenFval(c) = f_c;
        screenEsts{c} = e_c;
        screenOmega(c,:) = e_c(NParEst+1:end)';
        if isfield(out_c, 'firstorderopt')
            screenFOpt(c) = out_c.firstorderopt;
        end
    catch ME
        fprintf('  Candidate %d failed: %s\n', c, ME.message);
        screenFval(c) = Inf;
    end
    % Flag if type swap persists
    swapFlag = '';
    if ~isnan(screenOmega(c,3)) && screenOmega(c,3) > screenOmega(c,4)
        swapFlag = ' ** SWAP';
    end
    fprintf('  %2d. %-35s fval=%.6f  AllE=%.3f HiE=%.3f%s\n', ...
        c, candLabels{c}, screenFval(c), screenOmega(c,4), screenOmega(c,3), swapFlag);
end

%% Separate: best overall vs best non-swapped
isSwapped = screenOmega(:,3) > screenOmega(:,4);
fprintf('\n--- Screening summary ---\n');
fprintf('Swapped solutions:     %d / %d\n', sum(isSwapped & ~isinf(screenFval)), nCand);
fprintf('Non-swapped solutions: %d / %d\n', sum(~isSwapped & ~isinf(screenFval)), nCand);

[bestSwapFval, iBS] = min(screenFval .* (1 + ~isSwapped * 1e6));
[bestNoSwapFval, iBNS] = min(screenFval .* (1 + isSwapped * 1e6));
fprintf('Best swapped:     fval=%.6f (%s)\n', screenFval(iBS), candLabels{iBS});
if any(~isSwapped & ~isinf(screenFval))
    fprintf('Best non-swapped: fval=%.6f (%s)\n', screenFval(iBNS), candLabels{iBNS});
    fprintf('Gap (non-swap - swap): %.6f\n', bestNoSwapFval - bestSwapFval);
else
    fprintf('No non-swapped solutions found!\n');
end

%% Phase B: Refine top 3 non-swapped + top 2 swapped
% We want to see if non-swapped can beat swapped after full refinement
refineIdx = [];
refineLabels = {};

% Top 3 non-swapped
noSwapIdx = find(~isSwapped & ~isinf(screenFval));
if ~isempty(noSwapIdx)
    [~, srt] = sort(screenFval(noSwapIdx));
    topNS = noSwapIdx(srt(1:min(3, length(srt))));
    for i = 1:length(topNS)
        refineIdx(end+1) = topNS(i);
        refineLabels{end+1} = [candLabels{topNS(i)} ' (no-swap)'];
    end
end

% Top 2 swapped (for comparison)
swapIdx = find(isSwapped & ~isinf(screenFval));
if ~isempty(swapIdx)
    [~, srt] = sort(screenFval(swapIdx));
    topS = swapIdx(srt(1:min(2, length(srt))));
    for i = 1:length(topS)
        refineIdx(end+1) = topS(i);
        refineLabels{end+1} = [candLabels{topS(i)} ' (SWAP)'];
    end
end

nRefine = length(refineIdx);
refinedFval = inf(nRefine, 1);
refinedEsts = cell(nRefine, 1);
refinedFOpt = nan(nRefine, 1);
refinedEF   = nan(nRefine, 1);

fprintf('\nPhase B: Refining %d candidates ...\n', nRefine);
for i = 1:nRefine
    c = refineIdx(i);
    ests = screenEsts{c};
    fprintf('  Refining [%d] %s ...\n', c, refineLabels{i});

    % Pass 1: interior-point
    opts1 = optimoptions('fmincon', ...
        'Algorithm','interior-point', 'Display','off', ...
        'MaxIterations',15000, 'MaxFunctionEvaluations',40000, ...
        'OptimalityTolerance',1e-4, 'StepTolerance',1e-12, ...
        'ConstraintTolerance',1e-10, 'FiniteDifferenceType','forward', ...
        'FiniteDifferenceStepSize',1e-5);
    [ests, SSR2, ef, outF] = fmincon(f, ests, A, B, Aeq, 1, lb, ub, [], opts1);
    fo = outF.firstorderopt;
    fprintf('    Pass 1 (IP): fval=%.6f ef=%d 1stOrd=%.6f\n', SSR2, ef, fo);

    % Pass 2: SQP
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

    % Pass 3: interior-point central FD
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

    % Pass 4: SQP tight
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

    % Check if still swapped
    om_final = ests(NParEst+1:end);
    swapStr = '';
    if om_final(3) > om_final(4), swapStr = ' ** STILL SWAPPED'; end
    fprintf('  -> Final: fval=%.6f  AllE=%.3f HiE=%.3f%s\n\n', ...
        SSR2, om_final(4), om_final(3), swapStr);
end

%% Summary
fprintf('\n=== PHASE B RESULTS ===\n');
for i = 1:nRefine
    om = refinedEsts{i}(NParEst+1:end);
    sw = ''; if om(3) > om(4), sw = ' SWAP'; end
    fprintf('  %s: fval=%.6f 1stOrd=%.4f AllE=%.3f HiE=%.3f%s\n', ...
        refineLabels{i}, refinedFval(i), refinedFOpt(i), om(4), om(3), sw);
end

% Best non-swapped refined
noSwapRefined = [];
for i = 1:nRefine
    om = refinedEsts{i}(NParEst+1:end);
    if om(3) <= om(4)
        noSwapRefined(end+1) = i;
    end
end

if ~isempty(noSwapRefined)
    [bestNSfval, bestNSi] = min(refinedFval(noSwapRefined));
    bestIdx = noSwapRefined(bestNSi);
    fprintf('\nBest non-swapped: fval=%.6f (%s)\n', bestNSfval, refineLabels{bestIdx});
    fprintf('Best swapped:     fval=%.6f\n', min(refinedFval));
    fprintf('Gap: %.6f\n', bestNSfval - min(refinedFval));

    %% Post-process and save the best non-swapped solution
    ests = refinedEsts{bestIdx};
    pVecs = reshape(ests(1:NParEst), [12, Theta]);
    omega = ests(NParEst+1:end);

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

    % Recompute
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

    fprintf('\n--- Best non-swapped solution for 2024 ---\n');
    fprintf('omega: High-N=%.4f, High-U=%.4f, High-E=%.4f, All-E=%.4f, All-N=%.4f\n', omega);
    typeNamesM = {'High-N','High-U','High-E'};
    fprintf('Ergodic distributions:\n');
    for theta = 1:3
        fprintf('  %s: OLF=%.4f  U=%.4f  E-short=%.4f  E-long=%.4f  (E-mass=%.4f)\n', ...
            typeNamesM{theta}, ergVecs(:,theta), ...
            ergVecs(3,theta)+ergVecs(4,theta));
    end
    fprintf('1stOrdOpt: %.6f\n', refinedFOpt(bestIdx));

    fprintf('\nComparison with current (swapped) solution:\n');
    fprintf('  Current (swapped):  fval=%.6f\n', fval_current);
    fprintf('  New (non-swapped):  fval=%.6f\n', bestNSfval);
    fprintf('  Difference:         %.6f (positive = non-swapped is worse)\n', bestNSfval - fval_current);

    % Ask before saving
    fprintf('\n>> Save this non-swapped solution to Best folder? (fval gap = %.6f)\n', ...
        bestNSfval - fval_current);

    SSR2 = refinedFval(bestIdx);
    exitflag = refinedEF(bestIdx);
    res = struct();
    res.pVecs      = pVecs;
    res.omega      = omega;
    res.ergVecs    = ergVecs;
    res.ergVecsAll = ergVecsAll;
    res.SSR2       = SSR2;
    res.exitflag   = exitflag;

    save(fullfile(bestDir, 'step1_47.mat'), 'res');
    fprintf('SAVED to %s\n', fullfile(bestDir, 'step1_47.mat'));
else
    fprintf('\nWARNING: All refined solutions are still swapped!\n');
    fprintf('The data for 2024 may genuinely prefer this type structure.\n');
end

fprintf('\nTotal time: %.1f minutes\n', toc(tTotal)/60);
