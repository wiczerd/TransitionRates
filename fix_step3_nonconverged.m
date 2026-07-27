% fix_step3_nonconverged.m — Fix 4 non-converged years in Phase 3
%
% Years 2022-2024 barely iterated (33-38s, identical regPen=0.00122).
% Year 1987 ran 265s but didn't converge (regPen=0.00414).
%
% Strategy:
%   - 30 candidates per year (vs 9 in original run)
%   - Primary warm-starts: converged Phase 3 neighbor solutions (full 48 params)
%   - Secondary: Phase 1 warm-starts + All-E borrowed from converged neighbors
%   - Additional: wider random jitters, near-uniform, pure random
%   - Same 4-pass optimization structure

clc; clear;

%% ====================================================================
%  PATHS
%  ====================================================================
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(rewriteDir, '..', ...
    'Replication Files and ResultsTablesFigures Files 2 (1)');
phase1Dir  = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
step3Dir   = fullfile(rewriteDir, 'Step3_FixedOmega');

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

MaxEvals = 30000;
MaxIter  = 15000;
OptTols  = 0.001;

Nob        = 3^8;
Nstatepath = 4^8;
Nmonth     = 8;
Nact       = 3;
Nactpath   = Nact^Nmonth;

% Years to fix (indices into yearsList)
fixYears = [11, 45, 46, 47];  % 1987, 2022, 2023, 2024
fprintf('Fixing years: ');
for i = 1:numel(fixYears), fprintf('%d ', yearsList(fixYears(i))); end
fprintf('\n\n');

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
%  LOAD ALL PHASE 3 RESULTS (for warm-starts from converged years)
%  ====================================================================
step3_params = cell(numYears, 1);   % full 48-param vectors
step3_fvals  = nan(numYears, 1);
step3_ef     = nan(numYears, 1);
for yr = 1:numYears
    tmp = load(fullfile(step3Dir, sprintf('step3_%d.mat', yr)));
    step3_fvals(yr) = tmp.res.SSR2;
    step3_ef(yr)    = tmp.res.exitflag;
    % Reconstruct parameter vector from pVecs (already sorted)
    step3_params{yr} = tmp.res.pVecs(:);  % 48x1
end
converged = find(step3_ef == 1);
fprintf('Loaded %d Phase 3 results (%d converged)\n', numYears, numel(converged));

%% ====================================================================
%  LOAD PHASE 1 RESULTS (for additional warm-starts)
%  ====================================================================
phase1_pVecs = cell(numYears, 1);
for yr = 1:numYears
    tmp = load(fullfile(phase1Dir, sprintf('step1_%d.mat', yr)));
    phase1_pVecs{yr} = tmp.res.pVecs;  % 12x3
end
fprintf('Loaded Phase 1 warm-starts\n\n');

%% ====================================================================
%  FIX EACH YEAR
%  ====================================================================
epsIn = 1e-6;
allE_zeros = (Theta-1)*NParType + [2 5 7 8 9 11];  % structural zero positions

for fi = 1:numel(fixYears)
    year = fixYears(fi);
    tStart = tic;
    fprintf('=== Year %d (%d) ===\n', year, yearsList(year));
    fprintf('Original: fval=%.6f, ef=%d, regPen=%.2e\n', ...
        step3_fvals(year), step3_ef(year), NaN);

    Data_full = Data_full_TimeSeries(:, year);
    select = (1:Nob)';
    omeg1to5 = omega_fixed(year, 2:end)';

    f = @(x) ModelFitFunc_FixedOmega_AllE(x, select, groupNum, Theta, ...
        ThetaTot, Hstates, year, Data_full, omeg1to5);

    % ==================================================================
    % BUILD 30 CANDIDATES
    % ==================================================================
    candidates = {};
    cand_labels = {};

    % --- Group 1: Phase 3 converged neighbors (closest years first) ---
    % Sort all converged years by distance to this year
    dists = abs(converged - year);
    [~, sortI] = sort(dists);
    nearConv = converged(sortI);
    nNeighbors = min(10, numel(nearConv));
    for ni = 1:nNeighbors
        candidates{end+1} = step3_params{nearConv(ni)}; %#ok<AGROW>
        cand_labels{end+1} = sprintf('Phase3 yr%d', yearsList(nearConv(ni))); %#ok<AGROW>
    end

    % --- Group 2: Phase 1 warm-starts + All-E from nearest converged Phase 3 ---
    % Extract All-E params from the closest converged Phase 3 year
    allE_from_neighbor = step3_params{nearConv(1)}(37:48);  % last 12 params = All-E
    % This year's Phase 1
    p1_pv = phase1_pVecs{year};
    candidates{end+1} = [p1_pv(:); allE_from_neighbor];
    cand_labels{end+1} = sprintf('Phase1+AllE from yr%d', yearsList(nearConv(1)));

    % Phase 1 from ±1, ±2 neighbors
    neighbors = [year-1, year+1, year-2, year+2];
    neighbors = neighbors(neighbors >= 1 & neighbors <= numYears);
    for ni = 1:numel(neighbors)
        nb_pv = phase1_pVecs{neighbors(ni)};
        candidates{end+1} = [nb_pv(:); allE_from_neighbor]; %#ok<AGROW>
        cand_labels{end+1} = sprintf('Phase1 yr%d + AllE neighbor', yearsList(neighbors(ni))); %#ok<AGROW>
    end

    % --- Group 3: Phase 1 this year + generic All-E ---
    allE_init = [0.01; 0; 0.02; 0.01; 0; 0.01; 0; 0; 0; 0.01; 0; 0.98];
    candidates{end+1} = [p1_pv(:); allE_init];
    cand_labels{end+1} = 'Phase1 + generic AllE';

    % --- Group 4: Near-uniform ---
    unks_u = max(epsIn, min(1-epsIn, 0.25*0.9)) * ones(Nest, 1);
    unks_u(allE_zeros) = 0;
    candidates{end+1} = unks_u;
    cand_labels{end+1} = 'Near-uniform';

    % --- Group 5: Random candidates (diverse seeds) ---
    for rr = 1:5
        rng(100*year + rr, 'twister');
        Params_rand = [];
        for kk = 1:Theta
            p12 = rand; p13 = (1-p12)*rand; p14 = (1-p12-p13)*rand;
            p22 = rand; p23 = (1-p22)*rand; p24 = (1-p22-p23)*rand;
            p32 = rand; p33 = (1-p32)*rand; p34 = (1-p32-p33)*rand;
            p42 = rand; p43 = (1-p42)*rand; p44 = (1-p42-p43)*rand;
            if kk == Theta
                p13=0; p23=0; p32=0; p33=0; p34=0; p43=0;
                p14 = (1-p12)*rand; p24 = (1-p22)*rand; p44 = (1-p42)*rand;
            end
            Params_rand = [Params_rand; p12;p13;p14;p22;p23;p24;p32;p33;p34;p42;p43;p44]; %#ok<AGROW>
        end
        candidates{end+1} = Params_rand; %#ok<AGROW>
        cand_labels{end+1} = sprintf('Random seed %d', rr); %#ok<AGROW>
    end

    % --- Group 6: Jitters around best Phase 3 neighbors (wider: ±5%) ---
    for ji = 1:min(3, nNeighbors)
        for jj = 1:2
            rng(200*year + ji*10 + jj, 'twister');
            jitter = (rand(Nest, 1) - 0.5) * 0.10;  % ±5% absolute
            uJ = max(epsIn, min(1-epsIn, step3_params{nearConv(ji)} + jitter));
            uJ(allE_zeros) = 0;
            candidates{end+1} = uJ; %#ok<AGROW>
            cand_labels{end+1} = sprintf('Jitter(5%%) Phase3 yr%d #%d', yearsList(nearConv(ji)), jj); %#ok<AGROW>
        end
    end

    nCand = numel(candidates);
    fprintf('Built %d candidates\n', nCand);

    % ==================================================================
    % SCREEN CANDIDATES
    % ==================================================================
    cand_fvals = nan(nCand, 1);
    for c = 1:nCand
        cand_fvals(c) = f(candidates{c});
    end
    [sortedF, sortI] = sort(cand_fvals);

    fprintf('Top 5 candidates:\n');
    for k = 1:min(5, nCand)
        idx = sortI(k);
        fprintf('  %2d. fval=%.6f  %s\n', k, sortedF(k), cand_labels{idx});
    end

    % ==================================================================
    % OPTIMIZE TOP 5 CANDIDATES
    % ==================================================================
    nTop = min(5, nCand);
    topIdx = sortI(1:nTop);

    % Bounds
    lb = zeros(Nest, 1);
    ub = ones(Nest, 1);
    ub(allE_zeros) = 0;
    lb((Theta-1)*NParType + NParType) = 0.975;  % p44 for All-E

    % Row-sum constraints
    A_ineq = zeros(Theta*4, Nest);
    B_ineq = ones(Theta*4, 1);
    B_ineq((Theta-1)*4 + 3) = 0;  % All-E state 3 row = 0
    for index = 1:Theta*4
        A_ineq(index,:) = [zeros(1, 3*(index-1)), ones(1,3), ...
            zeros(1, Nest - 3 - 3*(index-1))];
    end

    best_ests = []; best_fval = inf; best_ef = 0; best_out = [];
    best_hessian = []; best_label = '';

    for ti = 1:nTop
        ci = topIdx(ti);
        unks = max(lb, min(ub, candidates{ci}));

        % Pass 1: interior-point, forward FD
        opts1 = optimoptions('fmincon', ...
            'Algorithm','interior-point', ...
            'Display','off', ...
            'MaxIterations',MaxIter, ...
            'MaxFunctionEvaluations',MaxEvals, ...
            'OptimalityTolerance',OptTols, ...
            'StepTolerance',1e-10, ...
            'ConstraintTolerance',1e-9, ...
            'FiniteDifferenceType','forward', ...
            'FiniteDifferenceStepSize',1e-5);
        [ests, SSR2, exitflag, outputF, ~, ~, hessian] = ...
            fmincon(f, unks, A_ineq, B_ineq, [], [], lb, ub, [], opts1);

        % Pass 2: SQP, central FD
        if exitflag == 2 || exitflag == 0 || ...
                (isfield(outputF,'firstorderopt') && outputF.firstorderopt > 1e-3)
            opts2 = optimoptions('fmincon', ...
                'Algorithm','sqp', ...
                'Display','off', ...
                'MaxIterations',MaxIter, ...
                'MaxFunctionEvaluations',MaxEvals, ...
                'OptimalityTolerance',OptTols, ...
                'StepTolerance',1e-10, ...
                'ConstraintTolerance',1e-9, ...
                'FiniteDifferenceType','central', ...
                'FiniteDifferenceStepSize',1e-4, ...
                'HessianApproximation','lbfgs');
            [e2, f2, ef2, out2, ~, ~, hes2] = ...
                fmincon(f, ests, A_ineq, B_ineq, [], [], lb, ub, [], opts2);
            if f2 < SSR2
                ests = min(ub, max(lb, e2));
                [SSR2, exitflag, outputF, hessian] = deal(f2, ef2, out2, hes2);
            end
        end

        % Pass 3: interior-point, central FD
        if exitflag == 2 || exitflag == 0 || ...
                (isfield(outputF,'firstorderopt') && outputF.firstorderopt > 1e-3)
            opts3 = optimoptions('fmincon', ...
                'Algorithm','interior-point', ...
                'Display','off', ...
                'MaxIterations',MaxIter, ...
                'MaxFunctionEvaluations',MaxEvals, ...
                'OptimalityTolerance',OptTols, ...
                'StepTolerance',1e-10, ...
                'ConstraintTolerance',1e-9, ...
                'FiniteDifferenceType','central', ...
                'FiniteDifferenceStepSize',1e-5);
            [e3, f3, ef3, out3, ~, ~, hes3] = ...
                fmincon(f, ests, A_ineq, B_ineq, [], [], lb, ub, [], opts3);
            if f3 < SSR2
                [ests, SSR2, exitflag, outputF, hessian] = deal(e3, f3, ef3, out3, hes3);
            end
        end

        % Pass 4: Polish
        ests_clean = min(ub, max(lb, ests));
        opts4 = optimoptions('fmincon', ...
            'Algorithm','interior-point', ...
            'Display','off', ...
            'MaxIterations',MaxIter, ...
            'MaxFunctionEvaluations',MaxEvals, ...
            'OptimalityTolerance',1e-6, ...
            'StepTolerance',1e-12, ...
            'ConstraintTolerance',1e-10, ...
            'FiniteDifferenceType','central', ...
            'FiniteDifferenceStepSize',1e-6);
        [e4, f4, ef4, out4, ~, ~, hes4] = ...
            fmincon(f, ests_clean, A_ineq, B_ineq, [], [], lb, ub, [], opts4);
        if f4 < SSR2 + 1e-8
            [ests, SSR2, exitflag, outputF, hessian] = deal(e4, f4, ef4, out4, hes4);
        end

        fo = 0;
        if isfield(outputF, 'firstorderopt'), fo = outputF.firstorderopt; end
        fprintf('  Candidate %d (%s): fval=%.6f, ef=%d, 1stOrd=%.2e\n', ...
            ti, cand_labels{ci}, SSR2, exitflag, fo);

        if SSR2 < best_fval
            best_ests = ests; best_fval = SSR2; best_ef = exitflag;
            best_out = outputF; best_hessian = hessian;
            best_label = cand_labels{ci};
        end
    end

    % ==================================================================
    % POST-PROCESSING (same as run_step3_fixedomega.m)
    % ==================================================================
    ests = best_ests;
    SSR2 = best_fval;
    exitflag = best_ef;
    outputF = best_out;
    hessian = best_hessian;

    omega = zeros(ThetaTot, 1);
    omega(1) = omeg1to5(2);
    omega(2) = 1 - omeg1to5(1) - omeg1to5(2) - omeg1to5(4) - omeg1to5(5);
    omega(3) = omeg1to5(4);
    omega(4) = omeg1to5(5);
    omega(5) = omeg1to5(1);

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

    % Type ordering: sort by E mass
    E_mass = ergVecs(3,:) + ergVecs(4,:);
    [~, I] = sort(E_mass);
    pVecs = pVecs(:, I);
    omega(1:Theta) = omega(I);

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

    % ==================================================================
    % DISPLAY ERGODIC DISTRIBUTIONS
    % ==================================================================
    typeNames = {'High-N','High-U','High-E','All-E'};
    stateNames = {'OLF','U','STJ','LTJ'};
    fprintf('\nErgodic distributions (after type ordering):\n');
    fprintf('  %-8s  %6s  %6s  %6s  %6s   omega\n', 'Type', stateNames{:});
    for theta = 1:Theta
        fprintf('  %-8s  %6.3f  %6.3f  %6.3f  %6.3f   %.4f\n', ...
            typeNames{theta}, ergVecs(:,theta)', omega(theta));
    end
    fprintf('  %-8s  %6.3f  %6.3f  %6.3f  %6.3f\n', 'All-N', 1, 0, 0, 0);
    fprintf('  %-8s  %6.3f  %6.3f  %6.3f  %6.3f   %.4f\n', ...
        'Aggregate', ergVecsAll', omega(ThetaTot));

    % ==================================================================
    % SAVE
    % ==================================================================
    fo = 0;
    if isfield(outputF, 'firstorderopt'), fo = outputF.firstorderopt; end

    res = struct();
    res.pVecs      = pVecs;
    res.omega      = omega;
    res.tms        = tms;
    res.ergVecs    = ergVecs;
    res.ergVecsAll = ergVecsAll;
    res.SSR2       = SSR2;
    res.exitflag   = exitflag;
    res.firstorder = fo;
    res.regPenalty = regPenalty;
    res.regTerms   = regTerms;
    res.regBinding = regBinding;

    % Compare with original
    improved = SSR2 < step3_fvals(year);
    fprintf('\nYear %d: fval %.6f -> %.6f (%s), ef=%d, regPen=%.2e, winner=%s\n', ...
        yearsList(year), step3_fvals(year), SSR2, ...
        ternary(improved, 'IMPROVED', 'NO CHANGE'), exitflag, regPenalty, best_label);

    if improved
        save(fullfile(step3Dir, sprintf('step3_%d.mat', year)), 'res');
        fprintf('Saved updated result.\n');
    else
        fprintf('Kept original result (no improvement).\n');
    end
    fprintf('Time: %.0fs\n\n', toc(tStart));
end

fprintf('Done.\n');

%% ====================================================================
function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
