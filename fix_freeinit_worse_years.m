% fix_freeinit_worse_years.m — Re-run the 8 years where free-init did worse
%
% These years found a worse local minimum with 50 params than the ergodic
% 41-param solution. Strategy: add more candidates including:
%   - The ergodic solution itself (pi = ergodic(tm)) as explicit start
%   - Free-init neighbors from already-completed years
%   - More jittered variants and random starts
%
% Saves improved results back to Step1_FreeOmega_FreeInitDist/

clc; clear all;
tTotal = tic;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(fileparts(rewriteDir), ...
               'Replication Files and ResultsTablesFigures Files 2 (1)');
bestDir    = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
fiDir      = fullfile(rewriteDir, 'Step1_FreeOmega_FreeInitDist');

addpath(origDir);
addpath(rewriteDir);

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

NParTM  = 12 * Theta;
NParPi  = 3 * Theta;
Nest    = NParTM + NParPi + ThetaTot;

% The 8 years where free-init was worse
badYears = [1986, 2003, 2008, 2009, 2011, 2012, 2022, 2024];
badIdx   = arrayfun(@(y) find(yearsList == y), badYears);
nBad     = length(badYears);

fprintf('=== Fixing %d years where free-init was worse ===\n', nBad);

%% Build Hstates
Hstates = zeros(Nstatepath, Nmonth);
path = 1;
for j1=1:4, for j2=1:4, for j3=1:4, for j4=1:4
for j5=1:4, for j6=1:4, for j7=1:4, for j8=1:4
    Hstates(path,:) = [j1 j2 j3 j4 j5 j6 j7 j8];
    path = path+1;
end, end, end, end, end, end, end, end

%% Load data
Data_full_TimeSeries = xlsread(fullfile(origDir, inputFileNames{groupNum}), 'J2:BD6562');

%% Load ergodic best-of-both results
fprintf('Loading ergodic best-of-both results ...\n');
bestPVecs  = cell(TotalYears, 1);
bestOmega  = cell(TotalYears, 1);
bestErgVec = cell(TotalYears, 1);
bestFval   = nan(TotalYears, 1);
for y = 1:TotalYears
    tmp = load(fullfile(bestDir, sprintf('step1_%d.mat', y)));
    bestPVecs{y}  = tmp.res.pVecs;
    bestOmega{y}  = tmp.res.omega;
    bestFval(y)   = tmp.res.SSR2;
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

%% Load free-init results (including neighbors that did well)
fprintf('Loading free-init results ...\n');
fiPVecs  = cell(TotalYears, 1);
fiOmega  = cell(TotalYears, 1);
fiInitD  = cell(TotalYears, 1);
fiFval   = nan(TotalYears, 1);
for y = 1:TotalYears
    tmp = load(fullfile(fiDir, sprintf('step1_%d.mat', y)));
    fiPVecs{y}  = tmp.res.pVecs;
    fiOmega{y}  = tmp.res.omega;
    fiInitD{y}  = tmp.res.initDists;
    fiFval(y)   = tmp.res.SSR2;
end

%% Helper: pack parameters
pack_unks = @(pVecs, piDist, omega) [pVecs(:); ...
    piDist(2,1); piDist(3,1); piDist(4,1); ...
    piDist(2,2); piDist(3,2); piDist(4,2); ...
    piDist(2,3); piDist(3,3); piDist(4,3); ...
    omega(:)];

%% Build constraints
lb = zeros(Nest, 1);
ub = ones(Nest, 1);
nIneq = Theta*4 + Theta;
A_ineq = zeros(nIneq, Nest);
B_ineq = ones(nIneq, 1);
for index = 1:Theta*4
    A_ineq(index, 3*(index-1)+1 : 3*index) = 1;
end
for theta = 1:Theta
    row = Theta*4 + theta;
    cols = NParTM + (theta-1)*3 + (1:3);
    A_ineq(row, cols) = 1;
end
Aeq = zeros(1, Nest);
Aeq(NParTM + NParPi + 1 : end) = 1;

%% Parallel pool
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local', min(nBad, 6));
end
fprintf('Parallel pool: %d workers\n', pool.NumWorkers);
pctRunOnAll addpath(origDir);
pctRunOnAll addpath(rewriteDir);

%% Run
jitterStd = 0.03;
Results = cell(nBad, 1);

parfor i = 1:nBad
    year = badIdx(i);
    Data_full = Data_full_TimeSeries(:, year);
    select = (1:3^8)';

    % =================================================================
    % BUILD CANDIDATES (~15 per year)
    % =================================================================
    candidates = {};
    nCand = 0;

    % 1. Ergodic solution with pi = ergodic(tm) (KEY ADDITION)
    nCand = nCand + 1;
    candidates{nCand} = pack_unks(bestPVecs{year}, bestErgVec{year}, bestOmega{year});

    % 2. Previous free-init result
    nCand = nCand + 1;
    candidates{nCand} = pack_unks(fiPVecs{year}, fiInitD{year}, fiOmega{year});

    % 3-4. Free-init neighbors t-1, t+1
    for nb = [-1, 1]
        ny = max(1, min(TotalYears, year + nb));
        nCand = nCand + 1;
        candidates{nCand} = pack_unks(fiPVecs{ny}, fiInitD{ny}, fiOmega{ny});
    end

    % 5-6. Ergodic neighbors t-1, t+1
    for nb = [-1, 1]
        ny = max(1, min(TotalYears, year + nb));
        nCand = nCand + 1;
        candidates{nCand} = pack_unks(bestPVecs{ny}, bestErgVec{ny}, bestOmega{ny});
    end

    % 7-10. Jittered ergodic solution (4 variants)
    for j = 1:4
        nCand = nCand + 1;
        rng(50000*year + j, 'twister');
        base = candidates{1};
        unks_j = base + jitterStd * randn(Nest, 1);
        unks_j = max(0, min(1, unks_j));
        om = unks_j(NParTM+NParPi+1:end);
        om = max(om, 1e-4);
        unks_j(NParTM+NParPi+1:end) = om / sum(om);
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

    % 11-14. Random starts
    for j = 1:4
        nCand = nCand + 1;
        rng(60000*year + j, 'twister');
        Params = [];
        for kk = 1:Theta
            p12 = rand; p13 = (1-p12)*rand; p14 = (1-p12-p13)*rand;
            p22 = rand; p23 = (1-p22)*rand; p24 = (1-p22-p23)*rand;
            p32 = rand; p33 = (1-p32)*rand; p34 = (1-p32-p33)*rand;
            p42 = rand; p43 = (1-p42)*rand; p44 = (1-p42-p43)*rand;
            Params = [Params; p12;p13;p14;p22;p23;p24;p32;p33;p34;p42;p43;p44];
        end
        piParams = [];
        for kk = 1:Theta
            r = -log(rand(4,1));
            r = r / sum(r);
            piParams = [piParams; r(2); r(3); r(4)];
        end
        omeg = (1/ThetaTot) * ones(ThetaTot, 1);
        candidates{nCand} = max(0, min(1, [Params; piParams; omeg]));
    end

    % =================================================================
    % OBJECTIVE
    % =================================================================
    f = @(x) ModelFitFuncTimeSeries_FreeInitDist( ...
        x, select, groupNum, Theta, Hstates, year, Data_full);

    % =================================================================
    % PHASE A: SCREEN
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
    opts1 = optimoptions('fmincon', ...
        'Algorithm','interior-point', 'Display','off', ...
        'MaxIterations',10000, 'MaxFunctionEvaluations',20000, ...
        'OptimalityTolerance',0.001, 'StepTolerance',1e-10, ...
        'ConstraintTolerance',1e-9, 'FiniteDifferenceType','forward', ...
        'FiniteDifferenceStepSize',1e-5);

    [ests, SSR2, exitflag, outputF] = fmincon(f, ests, A_ineq, B_ineq, Aeq, 1, lb, ub, [], opts1);

    if (exitflag==2 || (isfield(outputF,'firstorderopt') && outputF.firstorderopt>1e-3)) || exitflag==0
        opts2 = optimoptions(opts1, 'Algorithm','sqp', ...
            'FiniteDifferenceType','central', 'FiniteDifferenceStepSize',1e-4, ...
            'HessianApproximation','lbfgs', 'Display','off', 'ScaleProblem','obj-and-constr');
        [e2,f2,ef2,out2] = fmincon(f,ests,A_ineq,B_ineq,Aeq,1,lb,ub,[],opts2);
        if f2 < SSR2
            [ests,SSR2,exitflag,outputF] = deal(e2,f2,ef2,out2);
        end

        if exitflag==2 || exitflag==0
            opts3 = optimoptions('fmincon', ...
                'Algorithm','interior-point', 'Display','off', ...
                'MaxIterations',10000, 'MaxFunctionEvaluations',20000, ...
                'OptimalityTolerance',0.001, 'StepTolerance',1e-10, ...
                'ConstraintTolerance',1e-9, 'FiniteDifferenceType','central', ...
                'FiniteDifferenceStepSize',1e-5);
            [e3,f3,ef3,out3] = fmincon(f,ests,A_ineq,B_ineq,Aeq,1,lb,ub,[],opts3);
            if f3 < SSR2
                [ests,SSR2,exitflag,outputF] = deal(e3,f3,ef3,out3);
            end
        end
    end

    % =================================================================
    % POST-PROCESSING
    % =================================================================
    pVecs = [];
    omega = ests(end-(ThetaTot-1):end);
    for theta = 1:Theta
        pVec = ests((theta-1)*12+1:theta*12);
        pVecs = [pVecs, pVec];
    end

    initDists = zeros(4, Theta);
    for theta = 1:Theta
        piIdx = NParTM + (theta-1)*3 + (1:3);
        pi234 = ests(piIdx);
        initDists(:, theta) = [1-sum(pi234); pi234(:)];
    end

    % State ordering
    for theta = 1:Theta
        if pVecs(8,theta) > pVecs(12,theta)
            old_p = pVecs(:,theta);
            pVecs(:,theta) = [old_p(1);old_p(3);old_p(2); old_p(4);old_p(6);old_p(5); ...
                              old_p(10);old_p(12);old_p(11); old_p(7);old_p(9);old_p(8)];
            old_pi = initDists(:, theta);
            initDists(:, theta) = [old_pi(1); old_pi(2); old_pi(4); old_pi(3)];
        end
    end

    % Type ordering
    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs(:,theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        ergVecs(:,theta) = Ergodic(tm);
    end
    E_mass = ergVecs(3,:) + ergVecs(4,:);
    [~, I] = sort(E_mass);
    pVecs = pVecs(:, I);
    omega(1:Theta) = omega(I);
    initDists = initDists(:, I);

    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs(:,theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        ergVecs(:,theta) = Ergodic(tm);
    end

    U_mass = ergVecs(2, 1:2);
    [~, Isub] = sort(U_mass);
    pVecs(:,1:2) = pVecs(:, Isub);
    omega(1:2) = omega(Isub);
    initDists(:, 1:2) = initDists(:, Isub);

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

    fopt = NaN;
    if isfield(outputF, 'firstorderopt'), fopt = outputF.firstorderopt; end

    res = struct();
    res.pVecs      = pVecs;
    res.omega      = omega;
    res.initDists  = initDists;
    res.ergVecs    = ergVecs;
    res.ergVecsAll = ergVecsAll;
    res.SSR2       = SSR2;
    res.exitflag   = exitflag;
    res.bestCandID = bestC;
    res.screenFvals = screenFval(1:nCand);

    Results{i} = struct('year', year, 'yearCal', badYears(i), ...
        'fval', SSR2, 'exitflag', exitflag, 'fopt', fopt, ...
        'bestCand', bestC, 'res', res);

    fprintf('Year %d (%d): ef=%d, fval=%.4f, bestCand=%d, oldFI=%.4f, ergodic=%.4f\n', ...
        year, badYears(i), exitflag, SSR2, bestC, fiFval(year), bestFval(year));
end

%% Save improved results and summary
fprintf('\n=== Results ===\n');
fprintf('%-6s  %-10s  %-10s  %-10s  %-10s  %-6s\n', ...
    'Year', 'NewFval', 'OldFI', 'Ergodic', 'BestCand', 'Saved?');

for i = 1:nBad
    R = Results{i};
    year = R.year;
    oldFI = fiFval(year);
    erg   = bestFval(year);

    % Save if better than both old free-init AND ergodic
    newBest = R.fval < oldFI - 1e-8;
    saved = false;
    if newBest
        res = R.res;
        parsave_fix(fullfile(fiDir, sprintf('step1_%d.mat', year)), res);
        saved = true;
    end

    fprintf('%-6d  %-10.4f  %-10.4f  %-10.4f  %-10d  %s\n', ...
        R.yearCal, R.fval, oldFI, erg, R.bestCand, string(saved));
end

fprintf('\nTotal time: %.1f minutes\n', toc(tTotal)/60);

function parsave_fix(outFile, res)
    save(outFile, 'res');
end
