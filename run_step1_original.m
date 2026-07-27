% run_step1_original.m — Run Step 1 (free omegas) using the ORIGINAL code
%
% This is a thin wrapper around MultTypeEstTimeSeries_NewTrMBoot_20250919_V3.m
% that: (1) runs all 47 years in parallel, (2) saves to Step1_FreeOmega/,
% (3) exports to Excel and generates figures.

clc; clear all;
tTotal = tic;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
origDir = fullfile(fileparts(rewriteDir), ...
    'Replication Files and ResultsTablesFigures Files 2 (1)');
outDir = fullfile(rewriteDir, 'Step1_FreeOmega');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% Add original code directory to path (for objective fn + Ergodic.m)
addpath(origDir);

%% Settings (from original code)
TotalYears = 2024-1976+1-2;  % 47
inputFileNames = {'DataTimeSeriesYW.xlsx','DataTimeSeriesPW.xlsx', ...
                  'DataTimeSeriesYM.xlsx','DataTimeSeriesPM_1976_2023.xlsx'};
groupNum = 4;
Nrep = 1; rep = 1;
Nob = 3^8;
MaxEvals = 20000;
MaxIter  = 10000;
OptTols  = 0.001;
NumRespondents = [339039 1255294 344935 1164770];
s = 100;
Theta = 3;
ThetaTot = Theta + 2;
MaximumLikelihood = true;
Nstatepath = 4^8;
Nmonth = 8;
Nact = 3;
Nactpath = Nact^Nmonth;

yearsList = setdiff(1976:2024, [1977 1993]);
fprintf('=== Step 1: Free Omegas (original code) ===\n');
fprintf('Theta=%d mover + 2 polar = %d total, %d years\n', Theta, ThetaTot, TotalYears);

%% Build Hstates (original nested-loop order)
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
Iterations = nan(TotalYears, 1);
Seconds    = nan(TotalYears, 1);

%% Main estimation loop — using original code logic verbatim
fprintf('\n--- Starting estimation for %d years ---\n', TotalYears);

parfor year = 1:TotalYears
    tStart = tic;
    Data_full = Data_full_TimeSeries(:, year);

    % --- Original code: initial guess ---
    select = (1:3^8)';
    EstParams = ones(12,1);
    EstOmegas = ones(ThetaTot,1);

    Params = [];
    for kk = 1:Theta
        rng(s+kk)  % same seed as original
        p12 = rand(1,1);
        p13 = (1-p12).*rand(1,1);
        p14 = (1-p12-p13).*rand(1,1);
        p22 = rand(1,1);
        p23 = (1-p22).*rand(1,1);
        p24 = (1-p22-p23).*rand(1,1);
        p32 = rand(1,1);
        p33 = (1-p32).*rand(1,1);
        p34 = (1-p32-p33).*rand(1,1);
        p42 = rand(1,1);
        p43 = (1-p42).*rand(1,1);
        p44 = (1-p42-p43).*rand(1,1);
        temp0 = [p12 p13 p14 p22 p23 p24 p32 p33 p34 p42 p43 p44]';
        Params = [Params temp0]; %#ok<AGROW>
    end
    omeg1to5 = (1/ThetaTot)*ones(ThetaTot,1);

    % --- Build parameter vector ---
    ParMixVec = Params(:);
    EstParams3 = kron(ones(Theta,1), EstParams);
    NParEst = sum(EstParams3);
    EstParMix = [EstParams3; EstOmegas];
    ParMixVec = [ParMixVec; omeg1to5];

    unks = zeros(numel(EstParMix),1);
    loc = 1;
    for par = 1:numel(EstParMix)
        if EstParMix(par)
            unks(loc) = ParMixVec(par);
            unks(loc) = min(1, max(0, unks(loc)));
            loc = loc+1;
        end
    end

    % --- Constraints (original) ---
    Nest = size(unks,1);
    lb = zeros(Nest,1);
    ub = ones(Nest,1);
    Aeq = [zeros(1,NParEst), ones(1,ThetaTot)];
    A = zeros(Theta*4, NParEst+ThetaTot);
    B = ones(Theta*4, 1);
    for index = 1:Theta*4
        A(index,:) = [zeros(1,3*(index-1)) ones(1,3) zeros(1,NParEst+ThetaTot-3-3*(index-1))];
    end

    % --- Objective (original) ---
    f = @(x) ModelFitFuncTimeSeries_NewTrMBoot_20250919_V2(x, select, groupNum, Theta, Hstates, year, Data_full);

    % --- Optimizer pass 1: interior-point (original options) ---
    options = optimoptions('fmincon', ...
        'Algorithm','interior-point', 'Display','off', ...
        'MaxIterations',MaxIter, 'MaxFunctionEvaluations',MaxEvals, ...
        'OptimalityTolerance',OptTols, 'StepTolerance',1e-10, ...
        'ConstraintTolerance',1e-9, 'FiniteDifferenceType','forward', ...
        'FiniteDifferenceStepSize',1e-5);

    [ests, SSR2, exitflag, outputF, ~, ~, hessian] = ...
        fmincon(f, unks, A, B, Aeq, 1, lb, ub, [], options);

    % --- Pass 2: SQP refinement (original logic) ---
    if (exitflag==2 || (isfield(outputF,'firstorderopt') && outputF.firstorderopt>1e-3)) || exitflag==0
        opts2 = optimoptions(options, 'Algorithm','sqp', ...
            'FiniteDifferenceType','central', 'FiniteDifferenceStepSize',1e-4, ...
            'HessianApproximation','lbfgs', 'Display','off', 'ScaleProblem','obj-and-constr');
        [e2,f2,ef2,out2,~,~,hes2] = fmincon(f,ests,A,B,Aeq,1,lb,ub,[],opts2);
        if f2 < SSR2
            [ests,SSR2,exitflag,outputF,hessian] = deal(e2,f2,ef2,out2,hes2);
        end

        % --- Pass 3: interior-point central (original logic) ---
        if exitflag==2 || exitflag==0
            opts3 = optimoptions('fmincon', ...
                'Algorithm','interior-point', 'Display','off', ...
                'MaxIterations',MaxIter, 'MaxFunctionEvaluations',MaxEvals, ...
                'OptimalityTolerance',OptTols, 'StepTolerance',1e-10, ...
                'ConstraintTolerance',1e-9, 'FiniteDifferenceType','central', ...
                'FiniteDifferenceStepSize',1e-5);
            [e3,f3,ef3,out3,~,~,hes3] = fmincon(f,ests,A,B,Aeq,1,lb,ub,[],opts3);
            if f3 < SSR2
                [ests,SSR2,exitflag,outputF,hessian] = deal(e3,f3,ef3,out3,hes3);
            end
        end
    end

    % --- Post-estimation: state ordering + type relabeling (original) ---
    pVecs = [];
    omega = ests(end-(ThetaTot-1):end);
    for theta = 1:Theta
        pVec = ests((theta-1)*12+1:theta*12);
        pVecs = [pVecs, pVec]; %#ok<AGROW>
    end

    % State ordering: p33 <= p44
    for theta = 1:Theta
        if pVecs(8,theta) > pVecs(12,theta)
            old = pVecs(:,theta);
            pVecs(:,theta) = [old(1);old(3);old(2); old(4);old(6);old(5); ...
                              old(10);old(12);old(11); old(7);old(9);old(8)];
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

    % Type ordering: sort by ergodic E mass (states 3+4)
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

    % Sub-sort types 1 and 2 by U mass (original: ergVecs(2,1:2))
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

    % --- Pack results (original format) ---
    Valid = zeros(1, Theta);
    results = [pVecs; ...
        [SSR2, zeros(1,Theta-1)]; ...
        omega(1:Theta)'; ...
        [omega(ThetaTot-1), zeros(1,Theta-1)]; ...
        [omega(ThetaTot),   zeros(1,Theta-1)]; ...
        [OptTols,           zeros(1,Theta-1)]; ...
        [s,                 zeros(1,Theta-1)]; ...
        Valid; ...
        [exitflag,          zeros(1,Theta-1)]; ...
        ergVecs];

    % Save (original format + extras for figures)
    res = struct();
    res.results = results(:);
    res.pVecs = pVecs;
    res.omega = omega;
    res.ergVecs = ergVecs;
    res.ergVecsAll = ergVecsAll;
    res.SSR2 = SSR2;
    res.exitflag = exitflag;

    parsave_orig(fullfile(outDir, sprintf('step1_%d.mat', year)), res);

    % Track convergence
    ExitFlag(year) = exitflag;
    Fval(year) = SSR2;
    if isfield(outputF,'iterations'), Iterations(year) = outputF.iterations; end
    if isfield(outputF,'firstorderopt'), FirstOrder(year) = outputF.firstorderopt; end
    Seconds(year) = toc(tStart);

    fprintf('Year %d (%d): ef=%d, fval=%.4f, %.0fs\n', ...
        year, yearsList(year), exitflag, SSR2, Seconds(year));
end

%% Convergence summary
T = table(yearsList(:), ExitFlag, FirstOrder, Iterations, Fval, Seconds, ...
    'VariableNames', {'Year','ExitFlag','FirstOrderOpt','Iterations','Fval','Seconds'});
disp(T);
writetable(T, fullfile(outDir, 'convergence_step1.csv'));
save(fullfile(outDir, 'convergence_step1.mat'), 'T');

fprintf('\n=== Step 1 complete ===\n');
fprintf('ExitFlag=1: %d / %d years\n', sum(ExitFlag==1), TotalYears);
fprintf('Total time: %.1f minutes\n', toc(tTotal)/60);

%% Export to Excel and generate figures
fprintf('\nExporting to Excel...\n');
run(fullfile(rewriteDir, 'export_step1_orig.m'));
fprintf('Generating figures...\n');
run(fullfile(rewriteDir, 'figures_step1_orig.m'));

function parsave_orig(outFile, res)
    save(outFile, 'res');
end
