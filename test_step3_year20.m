% test_step3_year20.m — Test Phase 3 on year 20 (1997) only
%
% Quick single-year test to verify regularization + no ScaleProblem
% Before: fval=1.814 (ScaleProblem premature), 1.768 (weak reg), 1.525 (no reg, swapped)
% Goal:   fval near 1.52-1.55 with correctly ordered types

clc; clear;

rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(rewriteDir, '..', ...
    'Replication Files and ResultsTablesFigures Files 2 (1)');
phase1Dir  = fullfile(rewriteDir, 'Step1_FreeOmega_Best');

year = 20;  % 1997
numYears = 47;
yearsList = setdiff(1976:2024, [1977 1993]);
Theta     = 4;
ThetaTot  = 5;
NParType  = 12;
Nest      = Theta * NParType;
groupNum  = 4;

MaxEvals = 20000;
MaxIter  = 10000;
OptTols  = 0.001;
s        = 100;

Nob        = 3^8;
Nstatepath = 4^8;
Nmonth     = 8;
Nact       = 3;

% Build Hstates
Hstates = zeros(Nstatepath, Nmonth);
path = 1;
for j1=1:4, for j2=1:4, for j3=1:4, for j4=1:4
for j5=1:4, for j6=1:4, for j7=1:4, for j8=1:4
    Hstates(path,:) = [j1 j2 j3 j4 j5 j6 j7 j8];
    path = path + 1;
end, end, end, end, end, end, end, end

% Load data
fprintf('Loading data...\n');
Data_full_TimeSeries = xlsread( ...
    fullfile(origDir, 'DataTimeSeriesPM_1976_2023.xlsx'), 'J2:BD6562');
Data_full = Data_full_TimeSeries(:, year);

% Load omegas
load(fullfile(rewriteDir, 'omega_fixed.mat'), 'omega_fixed');
omeg1to5 = omega_fixed(year, 2:end)';
fprintf('Year %d (%d): omegas = [%.4f %.4f %.4f %.4f %.4f]\n', ...
    year, yearsList(year), omeg1to5);

% Load Phase 1 warm-starts
phase1_pVecs = cell(numYears, 1);
for yr = 1:numYears
    tmp = load(fullfile(phase1Dir, sprintf('step1_%d.mat', yr)));
    phase1_pVecs{yr} = tmp.res.pVecs;
end

select = (1:Nob)';

% Objective
f = @(x) ModelFitFunc_FixedOmega_AllE(x, select, groupNum, Theta, ...
    ThetaTot, Hstates, year, Data_full, omeg1to5);

% ======================================================================
% BUILD CANDIDATES (same as run_step3_fixedomega.m)
% ======================================================================
epsIn = 1e-6;
allE_zeros = (Theta-1)*NParType + [2 5 7 8 9 11];

candidates = {};

% Phase 1 warm-start for this year
p1_pv = phase1_pVecs{year};
allE_init = [0.01; 0; 0.02; 0.01; 0; 0.01; 0; 0; 0; 0.01; 0; 0.98];
candidates{end+1} = [p1_pv(:); allE_init];

% Neighbors
neighbors = [year-1, year+1, year-2, year+2];
neighbors = neighbors(neighbors >= 1 & neighbors <= numYears);
for ni = 1:numel(neighbors)
    nb_pv = phase1_pVecs{neighbors(ni)};
    candidates{end+1} = [nb_pv(:); allE_init];
end

% Uniform
unks_u = max(epsIn, min(1-epsIn, 0.25*0.9)) * ones(Nest, 1);
unks_u(allE_zeros) = 0;
candidates{end+1} = unks_u;

% Random
Params_rand = [];
for kk = 1:Theta
    rng(s + kk, 'twister');
    p12 = rand; p13 = (1-p12)*rand; p14 = (1-p12-p13)*rand;
    p22 = rand; p23 = (1-p22)*rand; p24 = (1-p22-p23)*rand;
    p32 = rand; p33 = (1-p32)*rand; p34 = (1-p32-p33)*rand;
    p42 = rand; p43 = (1-p42)*rand; p44 = (1-p42-p43)*rand;
    if kk == Theta
        p13=0; p23=0; p32=0; p33=0; p34=0; p43=0;
        p14 = (1-p12)*rand; p24 = (1-p22)*rand; p44 = (1-p42)*rand;
    end
    Params_rand = [Params_rand; p12;p13;p14;p22;p23;p24;p32;p33;p34;p42;p43;p44];
end
candidates{end+1} = Params_rand;

% Jitters
for j = 1:2
    rng(s + year*100 + j, 'twister');
    jitter = (rand(Nest, 1) - 0.5) * 0.04;
    uJ = max(epsIn, min(1-epsIn, candidates{1} + jitter));
    uJ(allE_zeros) = 0;
    candidates{end+1} = uJ;
end

% Screen
fprintf('\nScreening %d candidates:\n', numel(candidates));
bestx = candidates{1}; bestf = inf;
for c = 1:numel(candidates)
    ftry = f(candidates{c});
    fprintf('  Candidate %d: fval=%.4f\n', c, ftry);
    if ftry < bestf
        bestf = ftry; bestx = candidates{c};
    end
end
fprintf('Best candidate: fval=%.4f\n\n', bestf);
unks = bestx;

% ======================================================================
% BOUNDS AND CONSTRAINTS
% ======================================================================
lb = zeros(Nest, 1);
ub = ones(Nest, 1);
ub(allE_zeros) = 0;
lb((Theta-1)*NParType + NParType) = 0.975;

A_ineq = zeros(Theta*4, Nest);
B_ineq = ones(Theta*4, 1);
B_ineq((Theta-1)*4 + 3) = 0;
for index = 1:Theta*4
    A_ineq(index,:) = [zeros(1, 3*(index-1)), ones(1,3), zeros(1, Nest-3-3*(index-1))];
end

% ======================================================================
% PASS 1: interior-point, forward FD, NO ScaleProblem
% ======================================================================
fprintf('--- Pass 1: interior-point, forward FD ---\n');
opts1 = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','iter', ...
    'MaxIterations',MaxIter, ...
    'MaxFunctionEvaluations',MaxEvals, ...
    'OptimalityTolerance',OptTols, ...
    'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-9, ...
    'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',1e-5);

tic;
[ests, SSR2, exitflag, outputF, ~, ~, hessian] = ...
    fmincon(f, unks, A_ineq, B_ineq, [], [], lb, ub, [], opts1);
t1 = toc;
fprintf('Pass 1: ef=%d, fval=%.6f, 1stOrd=%.2e, %.1fs\n\n', ...
    exitflag, SSR2, outputF.firstorderopt, t1);

% ======================================================================
% PASS 2: SQP, central FD
% ======================================================================
if exitflag == 2 || exitflag == 0 || outputF.firstorderopt > 1e-3
    fprintf('--- Pass 2: SQP, central FD ---\n');
    opts2 = optimoptions('fmincon', ...
        'Algorithm','sqp', ...
        'Display','iter', ...
        'MaxIterations',MaxIter, ...
        'MaxFunctionEvaluations',MaxEvals, ...
        'OptimalityTolerance',OptTols, ...
        'StepTolerance',1e-10, ...
        'ConstraintTolerance',1e-9, ...
        'FiniteDifferenceType','central', ...
        'FiniteDifferenceStepSize',1e-4, ...
        'HessianApproximation','lbfgs');
    tic;
    [e2, f2, ef2, out2, ~, ~, hes2] = ...
        fmincon(f, ests, A_ineq, B_ineq, [], [], lb, ub, [], opts2);
    t2 = toc;
    fprintf('Pass 2: ef=%d, fval=%.6f, 1stOrd=%.2e, %.1fs\n', ...
        ef2, f2, out2.firstorderopt, t2);
    if f2 < SSR2
        ests = min(ub, max(lb, e2));
        [SSR2, exitflag, outputF, hessian] = deal(f2, ef2, out2, hes2);
        fprintf('  -> Accepted (improved)\n\n');
    else
        fprintf('  -> Rejected (no improvement)\n\n');
    end
end

% ======================================================================
% PASS 3: interior-point, central FD
% ======================================================================
if exitflag == 2 || exitflag == 0 || ...
        (isfield(outputF,'firstorderopt') && outputF.firstorderopt > 1e-3)
    fprintf('--- Pass 3: interior-point, central FD ---\n');
    opts3 = optimoptions('fmincon', ...
        'Algorithm','interior-point', ...
        'Display','iter', ...
        'MaxIterations',MaxIter, ...
        'MaxFunctionEvaluations',MaxEvals, ...
        'OptimalityTolerance',OptTols, ...
        'StepTolerance',1e-10, ...
        'ConstraintTolerance',1e-9, ...
        'FiniteDifferenceType','central', ...
        'FiniteDifferenceStepSize',1e-5);
    tic;
    [e3, f3, ef3, out3, ~, ~, hes3] = ...
        fmincon(f, ests, A_ineq, B_ineq, [], [], lb, ub, [], opts3);
    t3 = toc;
    fprintf('Pass 3: ef=%d, fval=%.6f, 1stOrd=%.2e, %.1fs\n', ...
        ef3, f3, out3.firstorderopt, t3);
    if f3 < SSR2
        [ests, SSR2, exitflag, outputF, hessian] = deal(e3, f3, ef3, out3, hes3);
        fprintf('  -> Accepted\n\n');
    else
        fprintf('  -> Rejected\n\n');
    end
end

% ======================================================================
% DIAGNOSTICS: ergodic distributions by type
% ======================================================================
fprintf('=== FINAL RESULT ===\n');
fprintf('fval = %.6f, exitflag = %d, 1stOrd = %.2e\n\n', ...
    SSR2, exitflag, outputF.firstorderopt);

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

fprintf('Type ergodic distributions (before state ordering):\n');
fprintf('%10s %8s %8s %8s %8s %8s\n', 'Type', 'OLF', 'U', 'STJ', 'LTJ', 'E-mass');
typeNames = {'High-N', 'High-U', 'High-E', 'All-E'};
for theta = 1:Theta
    p = pVecs(:, theta);
    tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
          1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
          1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
          1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
    if theta == Theta
        tm = [1-p(1)-0-p(3), p(1), 0, p(3); ...
              1-p(4)-0-p(6), p(4), 0, p(6); ...
              0, 0, 0, 0; ...
              1-p(10)-0-p(12), p(10), 0, p(12)];
        J = [1 2 4]; tm3 = tm(J,J);
        pi3 = Ergodic(tm3);
        ergVec = zeros(4,1); ergVec(J) = pi3;
    else
        ergVec = Ergodic(tm);
    end
    Emass = ergVec(3) + ergVec(4);
    fprintf('%10s %8.4f %8.4f %8.4f %8.4f %8.4f  omega=%.4f\n', ...
        typeNames{theta}, ergVec(1), ergVec(2), ergVec(3), ergVec(4), Emass, omega(theta));
end

fprintf('\nType ordering check:\n');
fprintf('  theta=1 (High-N): should dominate OLF (state 1)\n');
fprintf('  theta=2 (High-U): should dominate U (state 2)\n');
fprintf('  theta=3 (High-E): should dominate E (states 3+4)\n');
fprintf('  theta=4 (All-E):  structural zeros on state 3\n');
