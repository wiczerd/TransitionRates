% run_step1_freeinit.m — Step 1 with FREE initial distribution (fmincon multi-start)
%
% Same as rerun_all_multistart.m but uses ModelFitFuncTimeSeries_FreeInitDist
% which estimates pi_theta as free parameters instead of using ergodic(tm).
%
% Parameter layout (50 total):
%   1:36   = transition matrix params (12 per mover type x 3 types)
%   37:45  = initial distribution params (3 per mover type x 3 types)
%   46:50  = omegas (5 values, sum to 1)
%
% Starting points for pi params: ergodic distributions from best-of-both results.
%
% Output: Paper2_Rewrite/Step1_FreeOmega_FreeInitDist/

clc; clear all;
tTotal = tic;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(fileparts(rewriteDir), ...
               'Replication Files and ResultsTablesFigures Files 2 (1)');
bestDir    = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
outDir     = fullfile(rewriteDir, 'Step1_FreeOmega_FreeInitDist');
if ~exist(outDir, 'dir'), mkdir(outDir); end

addpath(origDir);
addpath(rewriteDir);  % for ModelFitFuncTimeSeries_FreeInitDist

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

% Parameter counts
NParTM  = 12 * Theta;          % 36 transition matrix params
NParPi  = 3 * Theta;           % 9 initial distribution params
Nest    = NParTM + NParPi + ThetaTot;  % 50 total

% Multi-start settings
nJitter    = 2;
jitterStd  = 0.03;
randomSeeds = [200, 300];
nRandom    = length(randomSeeds);

fprintf('=== Step 1 Free Init Dist (fmincon): %d years, %d params ===\n', TotalYears, Nest);

%% Build Hstates
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

%% Load best-of-both results (starting points)
fprintf('Loading best-of-both results from %s ...\n', bestDir);
bestPVecs  = cell(TotalYears, 1);
bestOmega  = cell(TotalYears, 1);
bestErgVec = cell(TotalYears, 1);  % ergodic dists to use as pi starting values
bestFval   = nan(TotalYears, 1);
for y = 1:TotalYears
    tmp = load(fullfile(bestDir, sprintf('step1_%d.mat', y)));
    bestPVecs{y}  = tmp.res.pVecs;
    bestOmega{y}  = tmp.res.omega;
    bestFval(y)    = tmp.res.SSR2;
    % Compute ergodic for each mover type -> starting pi
    ergV = zeros(4, Theta);
    for theta = 1:Theta
        p = tmp.res.pVecs(:, theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        ergV(:, theta) = Ergodic(tm);
    end
    bestErgVec{y} = ergV;
end
fprintf('  Best-of-both fval range: [%.4f, %.4f]\n', min(bestFval), max(bestFval));

%% Build constraints
% Lower/upper bounds
lb = zeros(Nest, 1);
ub = ones(Nest, 1);

% Inequality constraints: A*x <= B
% (a) TM row sums: p_{s,2} + p_{s,3} + p_{s,4} <= 1 (12 constraints)
% (b) Init dist sums: pi2 + pi3 + pi4 <= 1 (3 constraints)
nIneq = Theta*4 + Theta;
A_ineq = zeros(nIneq, Nest);
B_ineq = ones(nIneq, 1);

% TM row constraints (same as before)
for index = 1:Theta*4
    A_ineq(index, 3*(index-1)+1 : 3*index) = 1;
end

% Init dist constraints
for theta = 1:Theta
    row = Theta*4 + theta;
    cols = NParTM + (theta-1)*3 + (1:3);
    A_ineq(row, cols) = 1;
end

% Equality constraint: omegas sum to 1
Aeq = zeros(1, Nest);
Aeq(NParTM + NParPi + 1 : end) = 1;

%% Helper: pack parameters into unks vector
% unks = [pVecs(:); piParams; omega]
% piParams = [pi2_1, pi3_1, pi4_1, pi2_2, pi3_2, pi4_2, pi2_3, pi3_3, pi4_3]
pack_unks = @(pVecs, ergV, omega) [pVecs(:); ...
    ergV(2,1); ergV(3,1); ergV(4,1); ...
    ergV(2,2); ergV(3,2); ergV(4,2); ...
    ergV(2,3); ergV(3,3); ergV(4,3); ...
    omega(:)];

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
pctRunOnAll addpath(rewriteDir);
pctRunOnAll clear functions;

%% Pre-allocate
ExitFlag   = nan(TotalYears, 1);
Fval       = nan(TotalYears, 1);
FirstOrder = nan(TotalYears, 1);
BestCandID = nan(TotalYears, 1);
Seconds    = nan(TotalYears, 1);

%% ========== MAIN ESTIMATION LOOP ==========
fprintf('\n--- Starting free-init-dist estimation for %d years ---\n', TotalYears);

parfor year = 1:TotalYears
    tStart = tic;
    Data_full = Data_full_TimeSeries(:, year);
    select = (1:3^8)';

    % =================================================================
    % BUILD CANDIDATE INITIAL GUESSES
    % =================================================================
    candidates = {};
    nCand = 0;

    % 1. Best-of-both solution with its ergodic as pi
    nCand = nCand + 1;
    candidates{nCand} = pack_unks(bestPVecs{year}, bestErgVec{year}, bestOmega{year});

    % 2. Neighbor t-1
    if year > 1
        nCand = nCand + 1;
        candidates{nCand} = pack_unks(bestPVecs{year-1}, bestErgVec{year-1}, bestOmega{year-1});
    end

    % 3. Neighbor t+1
    if year < TotalYears
        nCand = nCand + 1;
        candidates{nCand} = pack_unks(bestPVecs{year+1}, bestErgVec{year+1}, bestOmega{year+1});
    end

    % 4-5. Jittered best solution
    for j = 1:nJitter
        nCand = nCand + 1;
        rng(10000*year + j, 'twister');
        base = candidates{1};  % best-of-both
        unks_j = base + jitterStd * randn(Nest, 1);
        unks_j = max(0, min(1, unks_j));
        % Renormalize omegas
        om = unks_j(NParTM+NParPi+1:end);
        om = max(om, 1e-4);
        unks_j(NParTM+NParPi+1:end) = om / sum(om);
        % Renormalize pi for each type
        for theta = 1:Theta
            piIdx = NParTM + (theta-1)*3 + (1:3);
            piVals = max(unks_j(piIdx), 1e-4);
            if sum(piVals) > 1 - 1e-4
                piVals = piVals * (1-1e-3) / sum(piVals);
            end
            unks_j(piIdx) = piVals;
        end
        candidates{nCand} = unks_j;
    end

    % 6-7. Random
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
            Params = [Params; p12;p13;p14;p22;p23;p24;p32;p33;p34;p42;p43;p44];
        end
        % Random initial dist on 4-simplex for each type
        piParams = [];
        for kk = 1:Theta
            rng(seed + 100 + kk, 'twister');
            r = -log(rand(4,1));  % Dirichlet(1,1,1,1)
            r = r / sum(r);
            piParams = [piParams; r(2); r(3); r(4)];
        end
        omeg = (1/ThetaTot) * ones(ThetaTot, 1);
        candidates{nCand} = max(0, min(1, [Params; piParams; omeg]));
    end

    % =================================================================
    % OBJECTIVE FUNCTION
    % =================================================================
    f = @(x) ModelFitFuncTimeSeries_FreeInitDist( ...
        x, select, groupNum, Theta, Hstates, year, Data_full);

    % =================================================================
    % PHASE A: QUICK SCREEN ALL CANDIDATES
    % =================================================================
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
            [e_c, f_c] = fmincon(f, candidates{c}, A_ineq, B_ineq, Aeq, 1, lb, ub, [], opts_screen);
            screenFval(c) = f_c;
            screenEsts{c} = e_c;
        catch
            screenFval(c) = Inf;
        end
    end

    [~, bestC] = min(screenFval);
    ests = screenEsts{bestC};
    SSR2 = screenFval(bestC);

    % =================================================================
    % PHASE B: FULL REFINEMENT (3-pass)
    % =================================================================

    % Pass 1: interior-point, forward FD
    opts1 = optimoptions('fmincon', ...
        'Algorithm','interior-point', 'Display','off', ...
        'MaxIterations',10000, 'MaxFunctionEvaluations',20000, ...
        'OptimalityTolerance',0.001, 'StepTolerance',1e-10, ...
        'ConstraintTolerance',1e-9, 'FiniteDifferenceType','forward', ...
        'FiniteDifferenceStepSize',1e-5);

    [ests, SSR2, exitflag, outputF, ~, ~, hessian] = ...
        fmincon(f, ests, A_ineq, B_ineq, Aeq, 1, lb, ub, [], opts1);

    % Pass 2: SQP with central FD
    if (exitflag==2 || (isfield(outputF,'firstorderopt') && outputF.firstorderopt>1e-3)) || exitflag==0
        opts2 = optimoptions(opts1, 'Algorithm','sqp', ...
            'FiniteDifferenceType','central', 'FiniteDifferenceStepSize',1e-4, ...
            'HessianApproximation','lbfgs', 'Display','off', 'ScaleProblem','obj-and-constr');
        [e2,f2,ef2,out2,~,~,hes2] = fmincon(f,ests,A_ineq,B_ineq,Aeq,1,lb,ub,[],opts2);
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
            [e3,f3,ef3,out3,~,~,hes3] = fmincon(f,ests,A_ineq,B_ineq,Aeq,1,lb,ub,[],opts3);
            if f3 < SSR2
                [ests,SSR2,exitflag,outputF,hessian] = deal(e3,f3,ef3,out3,hes3);
            end
        end
    end

    % =================================================================
    % POST-PROCESSING
    % =================================================================

    % Extract parameters
    pVecs = [];
    omega = ests(end-(ThetaTot-1):end);
    for theta = 1:Theta
        pVec = ests((theta-1)*12+1:theta*12);
        pVecs = [pVecs, pVec];
    end

    % Extract estimated initial distributions
    initDists = zeros(4, Theta);
    for theta = 1:Theta
        piIdx = NParTM + (theta-1)*3 + (1:3);
        pi234 = ests(piIdx);
        initDists(:, theta) = [1-sum(pi234); pi234(:)];
    end

    % State ordering: ensure p33 <= p44
    for theta = 1:Theta
        if pVecs(8,theta) > pVecs(12,theta)
            old_p = pVecs(:,theta);
            pVecs(:,theta) = [old_p(1);old_p(3);old_p(2); old_p(4);old_p(6);old_p(5); ...
                              old_p(10);old_p(12);old_p(11); old_p(7);old_p(9);old_p(8)];
            % Also swap init dist states 3 and 4
            old_pi = initDists(:, theta);
            initDists(:, theta) = [old_pi(1); old_pi(2); old_pi(4); old_pi(3)];
        end
    end

    % Ergodic distributions (for type labeling, not used in estimation)
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
    initDists = initDists(:, I);

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
    initDists(:, 1:2) = initDists(:, Isub);

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

    % =================================================================
    % SAVE
    % =================================================================
    res = struct();
    res.pVecs      = pVecs;
    res.omega      = omega;
    res.initDists  = initDists;   % NEW: estimated initial distributions
    res.ergVecs    = ergVecs;     % ergodic (for comparison)
    res.ergVecsAll = ergVecsAll;
    res.SSR2       = SSR2;
    res.exitflag   = exitflag;
    res.bestCandID   = bestC;
    res.screenFvals  = screenFval(1:nCand);

    parsave_fi(fullfile(outDir, sprintf('step1_%d.mat', year)), res);

    % Track convergence
    ExitFlag(year)   = exitflag;
    Fval(year)       = SSR2;
    BestCandID(year) = bestC;
    Seconds(year)    = toc(tStart);
    if isfield(outputF,'firstorderopt')
        FirstOrder(year) = outputF.firstorderopt;
    end

    fprintf('Year %d (%d): ef=%d, fval=%.4f, bestCand=%d, bestOfBoth=%.4f, %.0fs\n', ...
        year, yearsList(year), exitflag, SSR2, bestC, bestFval(year), Seconds(year));
end

%% ========== SUMMARY ==========
DeltaFval = Fval - bestFval;
T = table(yearsList(:), ExitFlag, FirstOrder, Fval, bestFval, DeltaFval, BestCandID, Seconds, ...
    'VariableNames', {'Year','ExitFlag','FirstOrderOpt','FreeInitFval','ErgodicFval','DeltaFval','BestCand','Seconds'});
disp(T);
writetable(T, fullfile(outDir, 'convergence_freeinit.csv'));
save(fullfile(outDir, 'convergence_freeinit.mat'), 'T');

nBetter = sum(Fval < bestFval - 1e-6);
nSame   = sum(abs(DeltaFval) <= 1e-6);
nWorse  = sum(Fval > bestFval + 1e-6);
fprintf('\n=== Free Init Dist (fmincon) Complete ===\n');
fprintf('Free init better: %d / %d years\n', nBetter, TotalYears);
fprintf('Same:             %d / %d years\n', nSame, TotalYears);
fprintf('Free init worse:  %d / %d years\n', nWorse, TotalYears);
fprintf('Total time: %.1f minutes\n', toc(tTotal)/60);

function parsave_fi(outFile, res)
    save(outFile, 'res');
end
