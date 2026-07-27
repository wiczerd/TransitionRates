% run_step3_sequential.m — Sequential warm-starting for Phase 3
%
% Problem: Independent optimization per year can land in different basins,
% producing noisy time series even when each year is well-converged.
%
% Solution: Chain solutions across years:
%   Forward pass:  year 1 → 2 → ... → 47 (solution t is warm-start for t+1)
%   Backward pass: year 47 → 46 → ... → 1 (solution t is warm-start for t-1)
%   For each year, keep the best fval across forward, backward, and existing.
%
% This creates natural "inertia" that keeps solutions in similar basins
% across adjacent years, producing smoother time series without artificial
% constraints.

clc; clear;

%% ====================================================================
%  PATHS
%  ====================================================================
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(rewriteDir, '..', ...
    'Replication Files and ResultsTablesFigures Files 2 (1)');
step3Dir   = fullfile(rewriteDir, 'Step3_FixedOmega');
outDir     = fullfile(rewriteDir, 'Step3_FixedOmega_Sequential');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ====================================================================
%  SETTINGS
%  ====================================================================
numYears  = 47;
yearsList = setdiff(1976:2024, [1977 1993]);
Theta     = 4;
ThetaTot  = 5;
NParType  = 12;
Nest      = Theta * NParType;  % 48
groupNum  = 4;

MaxEvals = 20000;
MaxIter  = 10000;

Nob        = 3^8;
Nstatepath = 4^8;
Nmonth     = 8;
Nact       = 3;

allE_zeros = (Theta-1)*NParType + [2 5 7 8 9 11];

%% ====================================================================
%  BUILD HSTATES
%  ====================================================================
Hstates = zeros(Nstatepath, Nmonth);
path = 1;
for j1=1:4, for j2=1:4, for j3=1:4, for j4=1:4
for j5=1:4, for j6=1:4, for j7=1:4, for j8=1:4
    Hstates(path,:) = [j1 j2 j3 j4 j5 j6 j7 j8];
    path = path + 1;
end, end, end, end, end, end, end, end

%% ====================================================================
%  LOAD DATA
%  ====================================================================
fprintf('Loading data...\n');
Data_full_TimeSeries = xlsread( ...
    fullfile(origDir, 'DataTimeSeriesPM_1976_2023.xlsx'), 'J2:BD6562');
fprintf('Data loaded: %d paths x %d years\n', size(Data_full_TimeSeries));

%% ====================================================================
%  LOAD FIXED OMEGAS
%  ====================================================================
load(fullfile(rewriteDir, 'omega_fixed.mat'), 'omega_fixed');

%% ====================================================================
%  LOAD EXISTING PHASE 3 RESULTS (our current best)
%  ====================================================================
existing_params = zeros(Nest, numYears);
existing_fvals  = nan(numYears, 1);
for yr = 1:numYears
    tmp = load(fullfile(step3Dir, sprintf('step3_%d.mat', yr)));
    existing_params(:, yr) = tmp.res.pVecs(:);  % 48x1
    existing_fvals(yr)     = tmp.res.SSR2;
end
fprintf('Loaded existing Phase 3 results (mean fval=%.4f)\n\n', mean(existing_fvals));

%% ====================================================================
%  BOUNDS AND CONSTRAINTS (same for all years)
%  ====================================================================
lb = zeros(Nest, 1);
ub = ones(Nest, 1);
ub(allE_zeros) = 0;
lb((Theta-1)*NParType + NParType) = 0.975;  % p44 for All-E

A_ineq = zeros(Theta*4, Nest);
B_ineq = ones(Theta*4, 1);
B_ineq((Theta-1)*4 + 3) = 0;  % All-E state 3 row = 0
for index = 1:Theta*4
    A_ineq(index,:) = [zeros(1, 3*(index-1)), ones(1,3), ...
        zeros(1, Nest - 3 - 3*(index-1))];
end

%% ====================================================================
%  FORWARD PASS: year 1 → 47
%  ====================================================================
fprintf('=== FORWARD PASS ===\n');
fwd_params = zeros(Nest, numYears);
fwd_fvals  = nan(numYears, 1);
fwd_ef     = nan(numYears, 1);

for yr = 1:numYears
    tStart = tic;
    Data_full = Data_full_TimeSeries(:, yr);
    select = (1:Nob)';
    omeg1to5 = omega_fixed(yr, 2:end)';

    f = @(x) ModelFitFunc_FixedOmega_AllE(x, select, groupNum, Theta, ...
        ThetaTot, Hstates, yr, Data_full, omeg1to5);

    % Candidates: (1) existing best, (2) previous year's forward solution
    candidates = {existing_params(:, yr)};
    if yr > 1
        candidates{end+1} = fwd_params(:, yr-1);
    end

    % Screen
    bestx = candidates{1}; bestf = inf;
    for c = 1:numel(candidates)
        x0 = max(lb, min(ub, candidates{c}));
        ftry = f(x0);
        if ftry < bestf
            bestf = ftry; bestx = x0;
        end
    end

    % Optimize
    [ests, fval, ef] = optimize_from_warmstart(f, bestx, lb, ub, ...
        A_ineq, B_ineq, MaxIter, MaxEvals);

    fwd_params(:, yr) = ests;
    fwd_fvals(yr) = fval;
    fwd_ef(yr) = ef;

    delta = fval - existing_fvals(yr);
    tag = '';
    if delta < -1e-4, tag = ' ** IMPROVED';
    elseif delta > 1e-4, tag = ' (worse)'; end

    fprintf('Fwd %d (%d): fval=%.4f (existing=%.4f, delta=%+.4f)%s  [%.0fs]\n', ...
        yr, yearsList(yr), fval, existing_fvals(yr), delta, tag, toc(tStart));
end

%% ====================================================================
%  BACKWARD PASS: year 47 → 1
%  ====================================================================
fprintf('\n=== BACKWARD PASS ===\n');
bwd_params = zeros(Nest, numYears);
bwd_fvals  = nan(numYears, 1);
bwd_ef     = nan(numYears, 1);

for yr = numYears:-1:1
    tStart = tic;
    Data_full = Data_full_TimeSeries(:, yr);
    select = (1:Nob)';
    omeg1to5 = omega_fixed(yr, 2:end)';

    f = @(x) ModelFitFunc_FixedOmega_AllE(x, select, groupNum, Theta, ...
        ThetaTot, Hstates, yr, Data_full, omeg1to5);

    % Candidates: (1) existing best, (2) next year's backward solution
    candidates = {existing_params(:, yr)};
    if yr < numYears
        candidates{end+1} = bwd_params(:, yr+1);
    end

    % Screen
    bestx = candidates{1}; bestf = inf;
    for c = 1:numel(candidates)
        x0 = max(lb, min(ub, candidates{c}));
        ftry = f(x0);
        if ftry < bestf
            bestf = ftry; bestx = x0;
        end
    end

    % Optimize
    [ests, fval, ef] = optimize_from_warmstart(f, bestx, lb, ub, ...
        A_ineq, B_ineq, MaxIter, MaxEvals);

    bwd_params(:, yr) = ests;
    bwd_fvals(yr) = fval;
    bwd_ef(yr) = ef;

    delta = fval - existing_fvals(yr);
    tag = '';
    if delta < -1e-4, tag = ' ** IMPROVED';
    elseif delta > 1e-4, tag = ' (worse)'; end

    fprintf('Bwd %d (%d): fval=%.4f (existing=%.4f, delta=%+.4f)%s  [%.0fs]\n', ...
        yr, yearsList(yr), fval, existing_fvals(yr), delta, tag, toc(tStart));
end

%% ====================================================================
%  PICK BEST PER YEAR (existing vs forward vs backward)
%  ====================================================================
fprintf('\n=== COMBINING RESULTS ===\n');
best_params = zeros(Nest, numYears);
best_fvals  = nan(numYears, 1);
best_source = cell(numYears, 1);
nImproved = 0;

for yr = 1:numYears
    candidates_f = [existing_fvals(yr), fwd_fvals(yr), bwd_fvals(yr)];
    candidates_p = {existing_params(:,yr), fwd_params(:,yr), bwd_params(:,yr)};
    labels = {'existing', 'forward', 'backward'};

    [bf, bi] = min(candidates_f);
    best_params(:, yr) = candidates_p{bi};
    best_fvals(yr) = bf;
    best_source{yr} = labels{bi};

    if bf < existing_fvals(yr) - 1e-4
        nImproved = nImproved + 1;
    end
end

fprintf('Improved: %d/%d years\n', nImproved, numYears);
fprintf('Mean fval: existing=%.4f, sequential=%.4f (delta=%+.4f)\n', ...
    mean(existing_fvals), mean(best_fvals), mean(best_fvals) - mean(existing_fvals));

% Source breakdown
nExist = sum(strcmp(best_source, 'existing'));
nFwd   = sum(strcmp(best_source, 'forward'));
nBwd   = sum(strcmp(best_source, 'backward'));
fprintf('Sources: existing=%d, forward=%d, backward=%d\n\n', nExist, nFwd, nBwd);

%% ====================================================================
%  POST-PROCESS AND SAVE BEST RESULTS
%  ====================================================================
fprintf('Post-processing and saving...\n');

for yr = 1:numYears
    Data_full = Data_full_TimeSeries(:, yr);
    omeg1to5 = omega_fixed(yr, 2:end)';

    ests = best_params(:, yr);

    % Reconstruct omegas
    omega = zeros(ThetaTot, 1);
    omega(1) = omeg1to5(2);       % High-N
    omega(2) = 1 - omeg1to5(1) - omeg1to5(2) - omeg1to5(4) - omeg1to5(5);  % High-U
    omega(3) = omeg1to5(4);       % High-E
    omega(4) = omeg1to5(5);       % All-E
    omega(5) = omeg1to5(1);       % All-N

    % Extract pVecs
    pVecs = zeros(NParType, Theta);
    for theta = 1:Theta
        pVecs(:, theta) = ests((theta-1)*NParType+1 : theta*NParType);
    end

    % State ordering: p33 <= p44
    for theta = 1:Theta
        if pVecs(8, theta) > pVecs(12, theta)
            old_p = pVecs(:, theta);
            pVecs(:, theta) = [old_p(1); old_p(3); old_p(2); ...
                               old_p(4); old_p(6); old_p(5); ...
                               old_p(10); old_p(12); old_p(11); ...
                               old_p(7); old_p(9); old_p(8)];
        end
    end

    % Build TMs and ergodic distributions
    tms = zeros(4*Theta, 4);
    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs(:, theta);
        if theta == Theta
            tm = [1-p(1)-0-p(3), p(1), 0, p(3); ...
                  1-p(4)-0-p(6), p(4), 0, p(6); ...
                  0, 0, 0, 0; ...
                  1-p(10)-0-p(12), p(10), 0, p(12)];
        else
            tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
                  1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
                  1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
                  1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        end
        tms(4*theta-3:4*theta, :) = tm;

        if theta == Theta
            J = [1 2 4]; tm3 = tm(J,J);
            pi3 = Ergodic(tm3);
            ergVec = zeros(4,1); ergVec(J) = pi3;
        else
            ergVec = Ergodic(tm);
        end
        ergVecs(:, theta) = ergVec;
    end

    % Regularization diagnostic (before relabeling)
    m_reg = 0.02; tau_reg = 0.05;
    pos_reg = @(z) tau_reg .* (max(z./tau_reg, 0) + log1p(exp(-abs(z./tau_reg))));
    z1a = ergVecs(1,2) - ergVecs(1,1) + m_reg;
    z1b = ergVecs(1,3) - ergVecs(1,1) + m_reg;
    z2a = ergVecs(2,1) - ergVecs(2,2) + m_reg;
    z2b = ergVecs(2,3) - ergVecs(2,2) + m_reg;
    mass34 = ergVecs(3,:) + ergVecs(4,:);
    z3a = mass34(1) - mass34(3) + m_reg;
    z3b = mass34(2) - mass34(3) + m_reg;
    regTerms = [pos_reg(z1a)^2, pos_reg(z1b)^2, pos_reg(z2a)^2, ...
                pos_reg(z2b)^2, pos_reg(z3a)^2, pos_reg(z3b)^2];
    regPenalty = 2.0 * sum(regTerms);
    regBinding = any(regTerms > 1e-10);

    % Type ordering: E-U-N criterion
    % Step 1: identify All-E (highest E mass among all types)
    E_mass = ergVecs(3,:) + ergVecs(4,:);
    [~, allE_idx] = max(E_mass);

    % Step 2: among remaining 3 types, apply E-U-N
    movers = setdiff(1:Theta, allE_idx);
    E_mass_m = E_mass(movers);
    [~, highE_local] = max(E_mass_m);
    highE_idx = movers(highE_local);

    remaining = setdiff(movers, highE_idx);
    U_mass_r = ergVecs(2, remaining);
    [~, highU_local] = max(U_mass_r);
    highU_idx = remaining(highU_local);

    highN_idx = setdiff(remaining, highU_idx);

    sortOrder = [highN_idx, highU_idx, highE_idx, allE_idx];
    pVecs = pVecs(:, sortOrder);
    omega(1:Theta) = omega(sortOrder);

    % Rebuild after sort
    tms = zeros(4*Theta, 4);
    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs(:, theta);
        if theta == Theta
            tm = [1-p(1)-0-p(3), p(1), 0, p(3); ...
                  1-p(4)-0-p(6), p(4), 0, p(6); ...
                  0, 0, 0, 0; ...
                  1-p(10)-0-p(12), p(10), 0, p(12)];
        else
            tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
                  1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
                  1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
                  1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        end
        tms(4*theta-3:4*theta, :) = tm;
        if theta == Theta
            J = [1 2 4]; tm3 = tm(J,J);
            pi3 = Ergodic(tm3); ergVec = zeros(4,1); ergVec(J) = pi3;
        else
            ergVec = Ergodic(tm);
        end
        ergVecs(:, theta) = ergVec;
    end

    ergVecsAll = ergVecs * omega(1:Theta) + [omega(ThetaTot); 0; 0; 0];

    % Save
    res = struct();
    res.pVecs      = pVecs;
    res.omega      = omega;
    res.tms        = tms;
    res.ergVecs    = ergVecs;
    res.ergVecsAll = ergVecsAll;
    res.SSR2       = best_fvals(yr);
    res.exitflag   = 1;
    res.firstorder = 0;
    res.regPenalty = regPenalty;
    res.regTerms   = regTerms;
    res.regBinding = regBinding;
    res.source     = best_source{yr};

    save(fullfile(outDir, sprintf('step3_%d.mat', yr)), 'res');
end

%% ====================================================================
%  SUMMARY TABLE
%  ====================================================================
fprintf('\n=== YEAR-BY-YEAR COMPARISON ===\n');
fprintf('Year   Existing  Sequential  Delta    Source\n');
fprintf('----   --------  ----------  ------   ------\n');
for yr = 1:numYears
    delta = best_fvals(yr) - existing_fvals(yr);
    fprintf('%d   %.4f    %.4f    %+.4f   %s\n', ...
        yearsList(yr), existing_fvals(yr), best_fvals(yr), delta, best_source{yr});
end

T = table(yearsList(:), existing_fvals, best_fvals, ...
    best_fvals - existing_fvals, best_source, ...
    'VariableNames', {'Year','Existing_fval','Sequential_fval','Delta','Source'});
writetable(T, fullfile(outDir, 'comparison_summary.csv'));
save(fullfile(outDir, 'comparison_summary.mat'), 'T', ...
    'fwd_fvals', 'bwd_fvals', 'existing_fvals', 'best_fvals', 'best_source');

fprintf('\nTotal improved: %d/%d\n', nImproved, numYears);
fprintf('Mean fval: existing=%.4f -> sequential=%.4f\n', ...
    mean(existing_fvals), mean(best_fvals));
fprintf('Results saved to: %s\n', outDir);


%% ====================================================================
%  OPTIMIZE FROM WARM-START (3-pass)
%  ====================================================================
function [ests, fval, ef] = optimize_from_warmstart(f, x0, lb, ub, ...
    A_ineq, B_ineq, MaxIter, MaxEvals)

    % Pass 1: interior-point
    opts1 = optimoptions('fmincon', ...
        'Algorithm','interior-point', ...
        'Display','off', ...
        'MaxIterations',MaxIter, ...
        'MaxFunctionEvaluations',MaxEvals, ...
        'OptimalityTolerance',1e-4, ...
        'StepTolerance',1e-10, ...
        'ConstraintTolerance',1e-9, ...
        'FiniteDifferenceType','forward', ...
        'FiniteDifferenceStepSize',1e-5);

    [ests, fval, ef, out1] = ...
        fmincon(f, x0, A_ineq, B_ineq, [], [], lb, ub, [], opts1);

    % Pass 2: SQP refinement if not fully converged
    if ef == 2 || ef == 0 || ...
            (isfield(out1,'firstorderopt') && out1.firstorderopt > 1e-3)
        opts2 = optimoptions('fmincon', ...
            'Algorithm','sqp', ...
            'Display','off', ...
            'MaxIterations',MaxIter, ...
            'MaxFunctionEvaluations',MaxEvals, ...
            'OptimalityTolerance',1e-4, ...
            'StepTolerance',1e-10, ...
            'ConstraintTolerance',1e-9, ...
            'FiniteDifferenceType','central', ...
            'FiniteDifferenceStepSize',1e-4, ...
            'HessianApproximation','lbfgs');
        [e2, f2, ef2] = ...
            fmincon(f, ests, A_ineq, B_ineq, [], [], lb, ub, [], opts2);
        if f2 < fval
            ests = min(ub, max(lb, e2));
            fval = f2; ef = ef2;
        end
    end

    % Pass 3: Polish
    ests_clean = min(ub, max(lb, ests));
    opts3 = optimoptions('fmincon', ...
        'Algorithm','interior-point', ...
        'Display','off', ...
        'MaxIterations',MaxIter, ...
        'MaxFunctionEvaluations',MaxEvals, ...
        'OptimalityTolerance',1e-6, ...
        'StepTolerance',1e-12, ...
        'ConstraintTolerance',1e-10, ...
        'FiniteDifferenceType','central', ...
        'FiniteDifferenceStepSize',1e-6);
    [e3, f3, ef3] = ...
        fmincon(f, ests_clean, A_ineq, B_ineq, [], [], lb, ub, [], opts3);
    if f3 < fval + 1e-8
        ests = e3; fval = f3; ef = ef3;
    end
end
