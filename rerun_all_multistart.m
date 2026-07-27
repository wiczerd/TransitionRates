% rerun_all_multistart.m — Re-run Step 1 for ALL 47 years with multi-start
%
% Strategy:
%   Phase A: Quick screen ~7 candidates per year (fast optimizer settings)
%   Phase B: Full 3-pass refinement on the best candidate (same as original)
%
% Candidates per year:
%   1. Old Step1 solution (self warm-start)
%   2. Neighbor t-1 solution
%   3. Neighbor t+1 solution
%   4-5. Jittered old solution (Gaussian noise sigma=0.03)
%   6-7. Random guesses (seeds 200, 300; original used 100)
%
% Output: Paper2_Rewrite/Step1_FreeOmega_Improving/step1_{1..47}.mat
% Old results in Step1_FreeOmega/ are NOT modified.

clc; clear all;
tTotal = tic;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(fileparts(rewriteDir), ...
               'Replication Files and ResultsTablesFigures Files 2 (1)');
oldDir     = fullfile(rewriteDir, 'Step1_FreeOmega');
outDir     = fullfile(rewriteDir, 'Step1_FreeOmega_Improving');
if ~exist(outDir, 'dir'), mkdir(outDir); end

addpath(origDir);

%% Settings (identical to original)
TotalYears = 47;
Theta    = 3;
ThetaTot = Theta + 2;
groupNum = 4;
inputFileNames = {'DataTimeSeriesYW.xlsx','DataTimeSeriesPW.xlsx', ...
                  'DataTimeSeriesYM.xlsx','DataTimeSeriesPM_1976_2023.xlsx'};
NumRespondents = [339039 1255294 344935 1164770];
Nstatepath = 4^8;
Nmonth     = 8;
yearsList  = setdiff(1976:2024, [1977 1993]);

% Multi-start settings
nJitter    = 2;          % jittered versions of old solution
jitterStd  = 0.03;       % std of Gaussian noise for jittering
randomSeeds = [200, 300]; % different from original seed=100
nRandom    = length(randomSeeds);

fprintf('=== Step 1 Multi-Start: %d years, ~7 candidates each ===\n', TotalYears);

%% Build Hstates (same nested-loop order as original)
Hstates = zeros(Nstatepath, Nmonth);
path = 1;
for j1=1:4, for j2=1:4, for j3=1:4, for j4=1:4
for j5=1:4, for j6=1:4, for j7=1:4, for j8=1:4
    Hstates(path,:) = [j1 j2 j3 j4 j5 j6 j7 j8];
    path = path+1;
end, end, end, end, end, end, end, end

%% Load data
fprintf('Loading data from %s ...\n', inputFileNames{groupNum});
Data_full_TimeSeries = xlsread(fullfile(origDir, inputFileNames{groupNum}), 'J2:BD6562');
fprintf('  Data: %d paths x %d years\n', size(Data_full_TimeSeries));

%% Load old Step1 results (to use as warm-starts)
fprintf('Loading old Step1 results from %s ...\n', oldDir);
oldUnks = cell(TotalYears, 1);
oldFval = nan(TotalYears, 1);
for y = 1:TotalYears
    tmp = load(fullfile(oldDir, sprintf('step1_%d.mat', y)));
    % Reconstruct unks vector: [pVecs(:); omega]
    oldUnks{y} = [tmp.res.pVecs(:); tmp.res.omega];
    oldFval(y) = tmp.res.SSR2;
end
fprintf('  Old fval range: [%.4f, %.4f]\n', min(oldFval), max(oldFval));

%% Build constraints (same for all years)
NParEst = 12 * Theta;          % 36 transition parameters
Nest    = NParEst + ThetaTot;  % 41 total parameters
lb  = zeros(Nest, 1);
ub  = ones(Nest, 1);
Aeq = [zeros(1, NParEst), ones(1, ThetaTot)];  % omegas sum to 1
A   = zeros(Theta*4, Nest);
B   = ones(Theta*4, 1);
for index = 1:Theta*4
    A(index,:) = [zeros(1,3*(index-1)) ones(1,3) zeros(1,Nest-3-3*(index-1))];
end

%% Parallel pool
pool = gcp('nocreate');
clust = parcluster('local');
maxW = clust.NumWorkers;
desiredW = min(TotalYears, maxW);
if isempty(pool)
    pool = parpool('local', desiredW);
elseif pool.NumWorkers ~= desiredW
    delete(pool); pool = parpool('local', desiredW);
end
fprintf('Parallel pool: %d workers\n', pool.NumWorkers);
pctRunOnAll addpath(origDir);
pctRunOnAll clear functions;
pctRunOnAll rng(0, 'twister');

%% Pre-allocate convergence tracking
ExitFlag   = nan(TotalYears, 1);
Fval       = nan(TotalYears, 1);
FirstOrder = nan(TotalYears, 1);
BestCandID = nan(TotalYears, 1);
Seconds    = nan(TotalYears, 1);
OldFvalVec = oldFval;  % for comparison in summary table

%% ========== MAIN ESTIMATION LOOP ==========
fprintf('\n--- Starting multi-start estimation for %d years ---\n', TotalYears);

parfor year = 1:TotalYears
    tStart = tic;
    Data_full = Data_full_TimeSeries(:, year);
    select = (1:3^8)';

    % =====================================================================
    % BUILD CANDIDATE INITIAL GUESSES
    % =====================================================================
    candidates = {};
    nCand = 0;

    % --- Candidate 1: old Step1 solution (self warm-start) ---
    nCand = nCand + 1;
    candidates{nCand} = oldUnks{year}; %#ok<PFBNS>

    % --- Candidate 2: neighbor year t-1 ---
    if year > 1
        nCand = nCand + 1;
        candidates{nCand} = oldUnks{year - 1};
    end

    % --- Candidate 3: neighbor year t+1 ---
    if year < TotalYears
        nCand = nCand + 1;
        candidates{nCand} = oldUnks{year + 1};
    end

    % --- Candidates 4-5: jittered old solution ---
    for j = 1:nJitter
        nCand = nCand + 1;
        rng(10000*year + j, 'twister');
        unks_j = oldUnks{year} + jitterStd * randn(Nest, 1);
        unks_j = max(0, min(1, unks_j));
        % Renormalize omegas to sum to 1
        om = unks_j(NParEst+1:end);
        om = max(om, 1e-4);  % keep strictly positive
        unks_j(NParEst+1:end) = om / sum(om);
        candidates{nCand} = unks_j;
    end

    % --- Candidates 6-7: random initial guesses (different seeds) ---
    for j = 1:nRandom
        nCand = nCand + 1;
        seed = randomSeeds(j);
        Params = [];
        for kk = 1:Theta
            rng(seed + kk, 'twister');
            p12 = rand; p13 = (1-p12)*rand; p14 = (1-p12-p13)*rand;
            p22 = rand; p23 = (1-p22)*rand; p24 = (1-p22-p23)*rand;
            p32 = rand; p33 = (1-p32)*rand; p34 = (1-p32-p33)*rand;
            p42 = rand; p43 = (1-p42)*rand; p44 = (1-p42-p43)*rand;
            Params = [Params; p12;p13;p14;p22;p23;p24;p32;p33;p34;p42;p43;p44]; %#ok<AGROW>
        end
        omeg = (1/ThetaTot) * ones(ThetaTot, 1);
        candidates{nCand} = max(0, min(1, [Params; omeg]));
    end

    % =====================================================================
    % OBJECTIVE FUNCTION (same as original)
    % =====================================================================
    f = @(x) ModelFitFuncTimeSeries_NewTrMBoot_20250919_V2( ...
        x, select, groupNum, Theta, Hstates, year, Data_full);

    % =====================================================================
    % PHASE A: QUICK SCREEN ALL CANDIDATES
    % Fast settings: fewer iterations, looser tolerance
    % =====================================================================
    opts_screen = optimoptions('fmincon', ...
        'Algorithm','interior-point', 'Display','off', ...
        'MaxIterations',3000, 'MaxFunctionEvaluations',8000, ...
        'OptimalityTolerance',0.01, 'StepTolerance',1e-8, ...
        'ConstraintTolerance',1e-9, 'FiniteDifferenceType','forward', ...
        'FiniteDifferenceStepSize',1e-5);

    screenFval = inf(nCand, 1);
    screenEsts = cell(nCand, 1);

    for c = 1:nCand
        try
            [e_c, f_c] = fmincon(f, candidates{c}, A, B, Aeq, 1, lb, ub, [], opts_screen);
            screenFval(c) = f_c;
            screenEsts{c} = e_c;
        catch
            screenFval(c) = Inf;
        end
    end

    % Pick the candidate with the lowest screening fval
    [~, bestC] = min(screenFval);
    ests = screenEsts{bestC};
    SSR2 = screenFval(bestC);

    % =====================================================================
    % PHASE B: FULL REFINEMENT ON BEST CANDIDATE (same 3-pass as original)
    % =====================================================================

    % Pass 1: interior-point, forward FD
    opts1 = optimoptions('fmincon', ...
        'Algorithm','interior-point', 'Display','off', ...
        'MaxIterations',10000, 'MaxFunctionEvaluations',20000, ...
        'OptimalityTolerance',0.001, 'StepTolerance',1e-10, ...
        'ConstraintTolerance',1e-9, 'FiniteDifferenceType','forward', ...
        'FiniteDifferenceStepSize',1e-5);

    [ests, SSR2, exitflag, outputF, ~, ~, hessian] = ...
        fmincon(f, ests, A, B, Aeq, 1, lb, ub, [], opts1);

    % Pass 2: SQP with central FD and L-BFGS Hessian
    if (exitflag==2 || (isfield(outputF,'firstorderopt') && outputF.firstorderopt>1e-3)) || exitflag==0
        opts2 = optimoptions(opts1, 'Algorithm','sqp', ...
            'FiniteDifferenceType','central', 'FiniteDifferenceStepSize',1e-4, ...
            'HessianApproximation','lbfgs', 'Display','off', 'ScaleProblem','obj-and-constr');
        [e2,f2,ef2,out2,~,~,hes2] = fmincon(f,ests,A,B,Aeq,1,lb,ub,[],opts2);
        if f2 < SSR2
            [ests,SSR2,exitflag,outputF,hessian] = deal(e2,f2,ef2,out2,hes2);
        end

        % Pass 3: interior-point with central FD
        if exitflag==2 || exitflag==0
            opts3 = optimoptions('fmincon', ...
                'Algorithm','interior-point', 'Display','off', ...
                'MaxIterations',10000, 'MaxFunctionEvaluations',20000, ...
                'OptimalityTolerance',0.001, 'StepTolerance',1e-10, ...
                'ConstraintTolerance',1e-9, 'FiniteDifferenceType','central', ...
                'FiniteDifferenceStepSize',1e-5);
            [e3,f3,ef3,out3,~,~,hes3] = fmincon(f,ests,A,B,Aeq,1,lb,ub,[],opts3);
            if f3 < SSR2
                [ests,SSR2,exitflag,outputF,hessian] = deal(e3,f3,ef3,out3,hes3);
            end
        end
    end

    % =====================================================================
    % POST-PROCESSING (identical to run_step1_original.m)
    % =====================================================================

    % --- Extract pVecs and omega from ests ---
    pVecs = [];
    omega = ests(end-(ThetaTot-1):end);
    for theta = 1:Theta
        pVec = ests((theta-1)*12+1:theta*12);
        pVecs = [pVecs, pVec]; %#ok<AGROW>
    end

    % --- State ordering: ensure p33 <= p44 (STJ vs LTJ) ---
    for theta = 1:Theta
        if pVecs(8,theta) > pVecs(12,theta)
            old_p = pVecs(:,theta);
            pVecs(:,theta) = [old_p(1);old_p(3);old_p(2); old_p(4);old_p(6);old_p(5); ...
                              old_p(10);old_p(12);old_p(11); old_p(7);old_p(9);old_p(8)];
        end
    end

    % --- Ergodic distributions ---
    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs(:,theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        ergVecs(:,theta) = Ergodic(tm);
    end

    % --- Type ordering: sort by ergodic E mass (states 3+4) ---
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

    % --- Sub-sort types 1 and 2 by U mass ---
    U_mass = ergVecs(2, 1:2);
    [~, Isub] = sort(U_mass);
    pVecs(:,1:2) = pVecs(:, Isub);
    omega(1:2) = omega(Isub);

    % --- Final ergodic ---
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

    % =====================================================================
    % PACK & SAVE (same format as run_step1_original.m + multi-start info)
    % =====================================================================
    Valid = zeros(1, Theta);
    results = [pVecs; ...
        [SSR2, zeros(1,Theta-1)]; ...
        omega(1:Theta)'; ...
        [omega(ThetaTot-1), zeros(1,Theta-1)]; ...
        [omega(ThetaTot),   zeros(1,Theta-1)]; ...
        [0.001,             zeros(1,Theta-1)]; ...
        [100,               zeros(1,Theta-1)]; ...
        Valid; ...
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
    res.bestCandID   = bestC;           % which candidate won
    res.screenFvals  = screenFval(1:nCand); % fval of all screened candidates

    parsave_ms(fullfile(outDir, sprintf('step1_%d.mat', year)), res);

    % Track convergence
    ExitFlag(year)   = exitflag;
    Fval(year)       = SSR2;
    BestCandID(year) = bestC;
    Seconds(year)    = toc(tStart);
    if isfield(outputF,'firstorderopt')
        FirstOrder(year) = outputF.firstorderopt;
    end

    fprintf('Year %d (%d): ef=%d, fval=%.4f, bestCand=%d, oldFval=%.4f, %.0fs\n', ...
        year, yearsList(year), exitflag, SSR2, bestC, OldFvalVec(year), Seconds(year));
end

%% ========== SUMMARY ==========
DeltaFval = Fval - OldFvalVec;
T = table(yearsList(:), ExitFlag, FirstOrder, Fval, OldFvalVec, DeltaFval, BestCandID, Seconds, ...
    'VariableNames', {'Year','ExitFlag','FirstOrderOpt','NewFval','OldFval','DeltaFval','BestCand','Seconds'});
disp(T);
writetable(T, fullfile(outDir, 'convergence_multistart.csv'));
save(fullfile(outDir, 'convergence_multistart.mat'), 'T');

nImproved = sum(Fval < OldFvalVec - 1e-6);
nConverged = sum(ExitFlag == 1);
fprintf('\n=== Multi-start complete ===\n');
fprintf('ExitFlag=1: %d / %d years\n', nConverged, TotalYears);
fprintf('Improved fval over old: %d / %d years (DeltaFval < 0)\n', nImproved, TotalYears);
fprintf('Total time: %.1f minutes\n', toc(tTotal)/60);

% Show the 6 previously-problematic years
badYears = [1984 1987 2003 2010 2011 2013];
fprintf('\n--- Previously problematic years ---\n');
for i = 1:length(badYears)
    idx = find(yearsList == badYears(i));
    fprintf('  %d (idx %d): old fval=%.4f -> new fval=%.4f, firstorderopt=%.6f, bestCand=%d\n', ...
        badYears(i), idx, OldFvalVec(idx), Fval(idx), FirstOrder(idx), BestCandID(idx));
end

function parsave_ms(outFile, res)
    save(outFile, 'res');
end
