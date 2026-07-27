% run_step3_fixedomega.m — Phase 3: Fixed-omega estimation
%
% Estimates year-by-year transition matrices with omegas fixed at their
% linear trend values (trend 3: linear excl 2020-2021).
%
% Phase 3 structure: Theta=4 non-polar types + 1 polar type (All-N)
%   theta=1: High-N  (regularization pushes to dominate state 1)
%   theta=2: High-U  (regularization pushes to dominate state 2)
%   theta=3: High-E  (regularization pushes to dominate states 3+4)
%   theta=4: All-E   (structural zeros on state 3, p44 >= 0.975)
%   type 5:  All-N   (polar, always OLF)
%
% 48 estimated parameters per year (4 types x 12 TM params each)
% Initial distribution: ergodic (will be relaxed in future runs)
%
% Key improvement: uses Phase 1 TM solutions as warm-starts for types 1-3.

clc; clear;

%% ====================================================================
%  PATHS
%  ====================================================================
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(rewriteDir, '..', ...
    'Replication Files and ResultsTablesFigures Files 2 (1)');
phase1Dir  = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
outDir     = fullfile(rewriteDir, 'Step3_FixedOmega');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ====================================================================
%  SETTINGS
%  ====================================================================
numYears  = 47;
yearsList = setdiff(1976:2024, [1977 1993]);
Theta     = 4;        % non-polar types (High-N, High-U, High-E, All-E)
ThetaTot  = 5;        % total types (+ All-N polar)
NParType  = 12;       % TM params per type
Nest      = Theta * NParType;  % 48 total params
groupNum  = 4;        % prime-age men

MaxEvals = 20000;
MaxIter  = 10000;
OptTols  = 0.001;
s        = 100;       % random seed

Nob        = 3^8;
Nstatepath = 4^8;
Nmonth     = 8;
Nact       = 3;
Nactpath   = Nact^Nmonth;

%% ====================================================================
%  BUILD HSTATES (65536 x 8)
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
fprintf('Loaded omega_fixed: %d years, trend 3 (linear excl 2020-2021)\n', ...
    size(omega_fixed, 1));

%% ====================================================================
%  LOAD PHASE 1 RESULTS FOR WARM-STARTS
%  ====================================================================
phase1_pVecs = cell(numYears, 1);
for yr = 1:numYears
    tmp = load(fullfile(phase1Dir, sprintf('step1_%d.mat', yr)));
    phase1_pVecs{yr} = tmp.res.pVecs;  % 12x3 [High-N, High-U, High-E]
end
fprintf('Loaded Phase 1 warm-starts for %d years\n', numYears);

%% ====================================================================
%  PARALLEL POOL
%  ====================================================================
pool = gcp('nocreate');
if isempty(pool)
    c = parcluster('Processes');
    nWorkers = min(numYears, c.NumWorkers);
    pool = parpool(nWorkers);
end
addAttachedFiles(pool, {'ModelFitFunc_FixedOmega_AllE.m', 'Ergodic.m'});
pctRunOnAll addpath(pwd)
pctRunOnAll clear functions
pctRunOnAll rng(0, 'twister')
fprintf('Parallel pool ready: %d workers\n\n', pool.NumWorkers);

%% ====================================================================
%  PREALLOCATE TRACKING
%  ====================================================================
ExitFlag   = nan(numYears, 1);
FirstOrder = nan(numYears, 1);
Fval       = nan(numYears, 1);
Seconds    = nan(numYears, 1);
RegPenalty = nan(numYears, 1);
RegBinding = nan(numYears, 1);

%% ====================================================================
%  MAIN ESTIMATION LOOP
%  ====================================================================
fprintf('Starting Phase 3 estimation: %d years, %d params/year\n', numYears, Nest);
fprintf('=========================================================\n');

parfor year = 1:numYears

    tStart = tic;
    Data_full = Data_full_TimeSeries(:, year);
    select = (1:Nob)';
    omeg1to5 = omega_fixed(year, 2:end)';

    % Objective function handle
    f = @(x) ModelFitFunc_FixedOmega_AllE(x, select, groupNum, Theta, ...
        ThetaTot, Hstates, year, Data_full, omeg1to5);

    % ==================================================================
    % BUILD CANDIDATES
    % ==================================================================
    epsIn = 1e-6;
    allE_zeros = (Theta-1)*NParType + [2 5 7 8 9 11];  % All-E structural zero positions

    candidates = {};

    % --- Candidate 1: Phase 1 warm-start for THIS year + sensible All-E ---
    p1_pv = phase1_pVecs{year};  % 12x3 from Phase 1 [High-N, High-U, High-E]
    allE_init = [0.01; 0; 0.02; 0.01; 0; 0.01; 0; 0; 0; 0.01; 0; 0.98];
    candidates{end+1} = [p1_pv(:); allE_init];

    % --- Candidates 2-5: Phase 1 warm-starts from NEIGHBORING years ---
    neighbors = [year-1, year+1, year-2, year+2];
    neighbors = neighbors(neighbors >= 1 & neighbors <= numYears);
    for ni = 1:numel(neighbors)
        nb_pv = phase1_pVecs{neighbors(ni)};
        candidates{end+1} = [nb_pv(:); allE_init]; %#ok<AGROW>
    end

    % --- Near-uniform interior (respects All-E zeros) ---
    unks_u = max(epsIn, min(1-epsIn, 0.25*0.9)) * ones(Nest, 1);
    unks_u(allE_zeros) = 0;
    candidates{end+1} = unks_u;

    % --- Random init (original method) ---
    Params_rand = [];
    for kk = 1:Theta
        rng(s + kk, 'twister');
        p12 = rand; p13 = (1-p12)*rand; p14 = (1-p12-p13)*rand;
        p22 = rand; p23 = (1-p22)*rand; p24 = (1-p22-p23)*rand;
        p32 = rand; p33 = (1-p32)*rand; p34 = (1-p32-p33)*rand;
        p42 = rand; p43 = (1-p42)*rand; p44 = (1-p42-p43)*rand;
        if kk == Theta  % All-E structural zeros
            p13=0; p23=0; p32=0; p33=0; p34=0; p43=0;
            p14 = (1-p12)*rand; p24 = (1-p22)*rand; p44 = (1-p42)*rand;
        end
        Params_rand = [Params_rand; p12;p13;p14;p22;p23;p24;p32;p33;p34;p42;p43;p44]; %#ok<AGROW>
    end
    candidates{end+1} = Params_rand;

    % --- Jitters around Phase 1 warm-start ---
    for j = 1:2
        rng(s + year*100 + j, 'twister');
        jitter = (rand(Nest, 1) - 0.5) * 0.04;  % +/-2% absolute
        uJ = max(epsIn, min(1-epsIn, candidates{1} + jitter));
        uJ(allE_zeros) = 0;
        candidates{end+1} = uJ; %#ok<AGROW>
    end

    % --- Screen candidates ---
    bestx = candidates{1}; bestf = inf;
    for c = 1:numel(candidates)
        ftry = f(candidates{c});
        if ftry < bestf
            bestf = ftry; bestx = candidates{c};
        end
    end
    unks = bestx;

    % ==================================================================
    % BOUNDS AND CONSTRAINTS
    % ==================================================================
    lb = zeros(Nest, 1);
    ub = ones(Nest, 1);
    ub(allE_zeros) = 0;                        % All-E structural zeros
    lb((Theta-1)*NParType + NParType) = 0.975;  % p44 for All-E near 1

    % Row-sum constraints: A*x <= B
    A_ineq = zeros(Theta*4, Nest);
    B_ineq = ones(Theta*4, 1);
    B_ineq((Theta-1)*4 + 3) = 0;  % All-E state 3 row = 0
    for index = 1:Theta*4
        A_ineq(index,:) = [zeros(1, 3*(index-1)), ones(1,3), ...
            zeros(1, Nest - 3 - 3*(index-1))];
    end

    % ==================================================================
    % PASS 1: INTERIOR-POINT, forward FD
    % ==================================================================
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

    % ==================================================================
    % PASS 2: SQP, central FD (if pass 1 not fully converged)
    % ==================================================================
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
            ests = min(ub, max(lb, e2));  % clip to bounds before next pass
            [SSR2, exitflag, outputF, hessian] = deal(f2, ef2, out2, hes2);
        end
    end

    % ==================================================================
    % PASS 3: INTERIOR-POINT, central FD (if still not converged)
    % ==================================================================
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

    % ==================================================================
    % PASS 4: POLISH — always run from best solution so far
    % ==================================================================
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
    if f4 < SSR2 + 1e-8  % accept if equal or better
        [ests, SSR2, exitflag, outputF, hessian] = deal(e4, f4, ef4, out4, hes4);
    end

    % ==================================================================
    % POST-PROCESSING
    % ==================================================================

    % --- Reconstruct omegas from fixed values ---
    omega = zeros(ThetaTot, 1);
    omega(1) = omeg1to5(2);  % High-N
    omega(2) = 1 - omeg1to5(1) - omeg1to5(2) - omeg1to5(4) - omeg1to5(5);  % High-U
    omega(3) = omeg1to5(4);  % High-E
    omega(4) = omeg1to5(5);  % All-E
    omega(5) = omeg1to5(1);  % All-N

    % --- Extract pVecs ---
    pVecs = zeros(NParType, Theta);
    for theta = 1:Theta
        pVecs(:, theta) = ests((theta-1)*NParType+1 : theta*NParType);
    end

    % --- State ordering: ensure p33 <= p44 (STJ vs LTJ) ---
    for theta = 1:Theta
        if pVecs(8, theta) > pVecs(12, theta)
            old_p = pVecs(:, theta);
            pVecs(:, theta) = [old_p(1); old_p(3); old_p(2); ...
                               old_p(4); old_p(6); old_p(5); ...
                               old_p(10); old_p(12); old_p(11); ...
                               old_p(7); old_p(9); old_p(8)];
        end
    end

    % --- Build transition matrices and ergodic distributions ---
    tms = zeros(4*Theta, 4);
    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs(:, theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];

        if theta == Theta  % All-E: enforce structural zeros
            tm = [1-p(1)-0-p(3), p(1), 0, p(3); ...
                  1-p(4)-0-p(6), p(4), 0, p(6); ...
                  0, 0, 0, 0; ...
                  1-p(10)-0-p(12), p(10), 0, p(12)];
        end
        tms(4*theta-3:4*theta, :) = tm;

        if theta == Theta
            J = [1 2 4];
            tm3 = tm(J, J);
            pi3 = Ergodic(tm3);
            ergVec = zeros(4, 1);
            ergVec(J) = pi3;
        else
            ergVec = Ergodic(tm);
        end
        ergVecs(:, theta) = ergVec;
    end

    % --- Regularization diagnostic (before any relabeling) ---
    % Compute individual penalty terms to check if regularization is binding
    m_reg = 0.02; tau_reg = 0.05;
    pos_reg = @(z) tau_reg .* (max(z./tau_reg, 0) + log1p(exp(-abs(z./tau_reg))));
    % theta=1 (High-N) should dominate OLF (state 1)
    z1a = ergVecs(1,2) - ergVecs(1,1) + m_reg;  % High-U OLF < High-N OLF
    z1b = ergVecs(1,3) - ergVecs(1,1) + m_reg;  % High-E OLF < High-N OLF
    % theta=2 (High-U) should dominate U (state 2)
    z2a = ergVecs(2,1) - ergVecs(2,2) + m_reg;  % High-N U < High-U U
    z2b = ergVecs(2,3) - ergVecs(2,2) + m_reg;  % High-E U < High-U U
    % theta=3 (High-E) should dominate E (states 3+4)
    mass34 = ergVecs(3,:) + ergVecs(4,:);
    z3a = mass34(1) - mass34(3) + m_reg;  % High-N E < High-E E
    z3b = mass34(2) - mass34(3) + m_reg;  % High-U E < High-E E
    regTerms = [pos_reg(z1a)^2, pos_reg(z1b)^2, pos_reg(z2a)^2, ...
                pos_reg(z2b)^2, pos_reg(z3a)^2, pos_reg(z3b)^2];
    regPenalty = 2.0 * sum(regTerms);  % lam * Nobs (Nobs≈1)
    regBinding = any(regTerms > 1e-10);

    % --- Type ordering: sort by ergodic E mass (states 3+4) ---
    E_mass = ergVecs(3,:) + ergVecs(4,:);
    [~, I] = sort(E_mass);
    pVecs = pVecs(:, I);
    omega(1:Theta) = omega(I);

    % Rebuild tms and ergVecs after sort
    tms = zeros(4*Theta, 4);
    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs(:, theta);
        if theta == Theta  % All-E (highest E mass → last after sort)
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

    % Aggregate ergodic (add All-N polar)
    ergVecsAll = ergVecs * omega(1:Theta) + [omega(ThetaTot); 0; 0; 0];

    % ==================================================================
    % SAVE RESULTS
    % ==================================================================
    res = struct();
    res.pVecs      = pVecs;       % 12 x 4 (ordered: High-N, High-U, High-E, All-E)
    res.omega      = omega;       % 5 x 1  (ordered: High-N, High-U, High-E, All-E, All-N)
    res.tms        = tms;         % 16 x 4 (stacked 4x4 transition matrices)
    res.ergVecs    = ergVecs;     % 4 x 4  (ergodic dist per type)
    res.ergVecsAll = ergVecsAll;  % 4 x 1  (aggregate ergodic)
    res.SSR2       = SSR2;        % scalar  (negative log-likelihood)
    res.exitflag   = exitflag;
    res.firstorder = 0;
    if isfield(outputF, 'firstorderopt')
        res.firstorder = outputF.firstorderopt;
    end
    res.regPenalty = regPenalty;   % total regularization penalty
    res.regTerms   = regTerms;    % 6 individual penalty terms
    res.regBinding = regBinding;  % true if any term > 1e-10

    parsave(fullfile(outDir, sprintf('step3_%d.mat', year)), res);

    % Track convergence
    ExitFlag(year)   = exitflag;
    Fval(year)       = SSR2;
    Seconds(year)    = toc(tStart);
    if isfield(outputF, 'firstorderopt')
        FirstOrder(year) = outputF.firstorderopt;
    end

    RegPenalty(year) = regPenalty;
    RegBinding(year) = regBinding;

    fprintf('Year %d (%d): ef=%d, fval=%.4f, 1stOrd=%.1e, regPen=%.2e, bind=%d, %.0fs\n', ...
        year, yearsList(year), exitflag, SSR2, FirstOrder(year), regPenalty, regBinding, Seconds(year));

end  % parfor year

%% ====================================================================
%  CONVERGENCE SUMMARY
%  ====================================================================
fprintf('\n=========================================================\n');
fprintf('Phase 3 estimation complete.\n\n');

T = table(yearsList(:), ExitFlag, FirstOrder, Fval, RegPenalty, RegBinding, Seconds, ...
    'VariableNames', {'Year','ExitFlag','FirstOrderOpt','Fval','RegPenalty','RegBinding','Seconds'});
disp(T);

writetable(T, fullfile(outDir, 'convergence_summary.csv'));
save(fullfile(outDir, 'convergence_summary.mat'), 'T');

nConverged = sum(ExitFlag == 1);
nRegBinding = sum(RegBinding == 1);
fprintf('Converged (ef=1): %d/%d\n', nConverged, numYears);
fprintf('Regularization binding: %d/%d years\n', nRegBinding, numYears);
fprintf('Total time: %.1f min\n', sum(Seconds)/60);
fprintf('Results saved to: %s\n', outDir);


%% ====================================================================
%  PARSAVE helper (for parfor)
%  ====================================================================
function parsave(fname, res)
    save(fname, 'res');
end
