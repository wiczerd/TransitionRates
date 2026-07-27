% run_step3_em_freeinit.m — Sequential EM with free initial distribution for Phase 3
%
% Same structure as run_step3_sequential.m (forward + backward + best-of-three)
% but at each year-step:
%   1. Run EM with free init dist (3 starts: existing sequential, prev year, random)
%   2. Refine with fmincon (free-init objective + regularization)
%   3. Pass solution to next year as warm-start
%
% Output: Step3_FixedOmega_FreeInit_EM/

clc; clear;
tTotal = tic;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(rewriteDir, '..', ...
    'Replication Files and ResultsTablesFigures Files 2 (1)');
seqDir     = fullfile(rewriteDir, 'Step3_FixedOmega_Sequential');
outDir     = fullfile(rewriteDir, 'Step3_FixedOmega_FreeInit_EM');
if ~exist(outDir, 'dir'), mkdir(outDir); end

addpath(origDir);
addpath(rewriteDir);

%% Settings
numYears  = 47;
yearsList = setdiff(1976:2024, [1977 1993]);
Theta     = 4;
ThetaTot  = 5;
NParType  = 12;
NParTM    = Theta * NParType;  % 48
NParPi    = 3 * Theta;          % 12
Nest      = NParTM + NParPi;    % 60
groupNum  = 4;

MaxEvals = 20000;
MaxIter  = 10000;

Nstatepath = 4^8;
Nmonth     = 8;
Nact       = 3;
Nactpath   = Nact^Nmonth;

% EM settings
maxEMiter   = 500;
emTol       = 1e-8;
pseudoCount = 1e-10;

fprintf('=== Phase 3 Sequential EM Free Init: %d years ===\n', numYears);

%% Build Hstates
Hstates = zeros(Nstatepath, Nmonth);
path = 1;
for j1=1:4, for j2=1:4, for j3=1:4, for j4=1:4
for j5=1:4, for j6=1:4, for j7=1:4, for j8=1:4
    Hstates(path,:) = [j1 j2 j3 j4 j5 j6 j7 j8];
    path = path+1;
end, end, end, end, end, end, end, end

%% Precompute mappings
map = [3, 2, 1, 1];
place_mults = Nact.^(7:-1:0)';
act_idx = zeros(Nstatepath, 1);
for p = 1:Nstatepath
    A = map(Hstates(p, :));
    act_idx(p) = (A - 1) * place_mults + 1;
end

h4h5_lin = sub2ind([4,4], Hstates(:,4), Hstates(:,5));
h4h5_groups = cell(16, 1);
for idx = 1:16
    h4h5_groups{idx} = find(h4h5_lin == idx);
end

h1_groups = cell(4, 1);
for s = 1:4
    h1_groups{s} = find(Hstates(:,1) == s);
end

direct_pairs = [1,2; 2,3; 3,4; 5,6; 6,7; 7,8];

%% Polar type
All_N = zeros(Nactpath, 1); All_N(end) = 1;

%% Load data
fprintf('Loading data ...\n');
Data_full_TimeSeries = xlsread(fullfile(origDir, ...
    'DataTimeSeriesPM_1976_2023.xlsx'), 'J2:BD6562');
fprintf('Data loaded: %d paths x %d years\n', size(Data_full_TimeSeries));

%% Load fixed omegas
load(fullfile(rewriteDir, 'omega_fixed.mat'), 'omega_fixed');

%% Load existing sequential results
fprintf('Loading sequential results ...\n');
existing_pVecs = cell(numYears, 1);
existing_fvals = nan(numYears, 1);
existing_ergV  = cell(numYears, 1);

for y = 1:numYears
    tmp = load(fullfile(seqDir, sprintf('step3_%d.mat', y)));
    existing_pVecs{y} = tmp.res.pVecs;
    existing_fvals(y) = tmp.res.SSR2;
    existing_ergV{y}  = tmp.res.ergVecs;
end
fprintf('Sequential results loaded (mean fval=%.4f)\n\n', mean(existing_fvals));

%% Bounds and constraints
lb = zeros(Nest, 1);
ub = ones(Nest, 1);
% All-E TM structural zeros
allE_zeros_tm = (Theta-1)*NParType + [2 5 7 8 9 11];
ub(allE_zeros_tm) = 0;
% All-E p44 lower bound
lb((Theta-1)*NParType + NParType) = 0.975;
% All-E pi(3) = 0
allE_pi3_idx = NParTM + (Theta-1)*3 + 2;
ub(allE_pi3_idx) = 0;

nTMConstr = Theta * 4;
nPiConstr = Theta;
nIneq = nTMConstr + nPiConstr;
A_ineq = zeros(nIneq, Nest);
B_ineq = ones(nIneq, 1);
B_ineq((Theta-1)*4 + 3) = 0;

for index = 1:nTMConstr
    A_ineq(index, 3*(index-1)+1 : 3*index) = 1;
end
for theta = 1:Theta
    row = nTMConstr + theta;
    cols = NParTM + (theta-1)*3 + (1:3);
    A_ineq(row, cols) = 1;
end

%% Pack shared constants into struct (for local functions)
cfg = struct();
cfg.Theta = Theta;
cfg.ThetaTot = ThetaTot;
cfg.Nstatepath = Nstatepath;
cfg.Nactpath = Nactpath;
cfg.NParTM = NParTM;
cfg.NParType = NParType;
cfg.groupNum = groupNum;
cfg.maxEMiter = maxEMiter;
cfg.emTol = emTol;
cfg.pseudoCount = pseudoCount;
cfg.MaxIter = MaxIter;
cfg.MaxEvals = MaxEvals;
cfg.Hstates = Hstates;
cfg.act_idx = act_idx;
cfg.h4h5_groups = h4h5_groups;
cfg.h1_groups = h1_groups;
cfg.direct_pairs = direct_pairs;
cfg.All_N = All_N;
cfg.lb = lb;
cfg.ub = ub;
cfg.A_ineq = A_ineq;
cfg.B_ineq = B_ineq;

%% ====================================================================
%  FORWARD PASS: year 1 -> 47
%  ====================================================================
fprintf('=== FORWARD PASS ===\n');
fwd_pVecs = cell(numYears, 1);
fwd_pi    = cell(numYears, 1);
fwd_fvals = nan(numYears, 1);

for yr = 1:numYears
    tStart = tic;
    M_data = Data_full_TimeSeries(:, yr);
    omeg1to5 = omega_fixed(yr, 2:end)';

    omega = zeros(ThetaTot, 1);
    omega(1) = omeg1to5(2);
    omega(2) = 1 - omeg1to5(1) - omeg1to5(2) - omeg1to5(4) - omeg1to5(5);
    omega(3) = omeg1to5(4);
    omega(4) = omeg1to5(5);
    omega(5) = omeg1to5(1);

    % Candidates: (1) existing sequential, (2) prev year fwd, (3) random
    cand_pV = {existing_pVecs{yr}};
    cand_pi = {existing_ergV{yr}};

    if yr > 1
        cand_pV{end+1} = fwd_pVecs{yr-1};
        cand_pi{end+1} = fwd_pi{yr-1};
    end

    rng(90000*yr + 1, 'twister');
    [rPV, rPi] = make_random_start(cfg);
    cand_pV{end+1} = rPV;
    cand_pi{end+1} = rPi;

    [pV, pi_v, fv] = run_year_em_fmincon(yr, M_data, omeg1to5, omega, cand_pV, cand_pi, cfg);
    fwd_pVecs{yr} = pV;
    fwd_pi{yr}    = pi_v;
    fwd_fvals(yr) = fv;

    delta = fv - existing_fvals(yr);
    tag = '';
    if delta < -1e-4, tag = ' ** IMPROVED';
    elseif delta > 1e-4, tag = ' (worse)'; end

    fprintf('Fwd %d (%d): fval=%.4f (seq=%.4f, delta=%+.4f)%s  [%.0fs]\n', ...
        yr, yearsList(yr), fv, existing_fvals(yr), delta, tag, toc(tStart));
end

%% ====================================================================
%  BACKWARD PASS: year 47 -> 1
%  ====================================================================
fprintf('\n=== BACKWARD PASS ===\n');
bwd_pVecs = cell(numYears, 1);
bwd_pi    = cell(numYears, 1);
bwd_fvals = nan(numYears, 1);

for yr = numYears:-1:1
    tStart = tic;
    M_data = Data_full_TimeSeries(:, yr);
    omeg1to5 = omega_fixed(yr, 2:end)';

    omega = zeros(ThetaTot, 1);
    omega(1) = omeg1to5(2);
    omega(2) = 1 - omeg1to5(1) - omeg1to5(2) - omeg1to5(4) - omeg1to5(5);
    omega(3) = omeg1to5(4);
    omega(4) = omeg1to5(5);
    omega(5) = omeg1to5(1);

    cand_pV = {existing_pVecs{yr}};
    cand_pi = {existing_ergV{yr}};

    if yr < numYears
        cand_pV{end+1} = bwd_pVecs{yr+1};
        cand_pi{end+1} = bwd_pi{yr+1};
    end

    rng(90000*yr + 2, 'twister');
    [rPV, rPi] = make_random_start(cfg);
    cand_pV{end+1} = rPV;
    cand_pi{end+1} = rPi;

    [pV, pi_v, fv] = run_year_em_fmincon(yr, M_data, omeg1to5, omega, cand_pV, cand_pi, cfg);
    bwd_pVecs{yr} = pV;
    bwd_pi{yr}    = pi_v;
    bwd_fvals(yr) = fv;

    delta = fv - existing_fvals(yr);
    tag = '';
    if delta < -1e-4, tag = ' ** IMPROVED';
    elseif delta > 1e-4, tag = ' (worse)'; end

    fprintf('Bwd %d (%d): fval=%.4f (seq=%.4f, delta=%+.4f)%s  [%.0fs]\n', ...
        yr, yearsList(yr), fv, existing_fvals(yr), delta, tag, toc(tStart));
end

%% ====================================================================
%  PICK BEST PER YEAR
%  ====================================================================
fprintf('\n=== COMBINING RESULTS ===\n');
best_pVecs = cell(numYears, 1);
best_pi    = cell(numYears, 1);
best_fvals = nan(numYears, 1);
best_source = cell(numYears, 1);
nImproved = 0;

for yr = 1:numYears
    candidates_f = [existing_fvals(yr), fwd_fvals(yr), bwd_fvals(yr)];
    candidates_pV = {existing_pVecs{yr}, fwd_pVecs{yr}, bwd_pVecs{yr}};
    candidates_pi = {existing_ergV{yr}, fwd_pi{yr}, bwd_pi{yr}};
    labels = {'existing', 'forward', 'backward'};

    [bf, bi] = min(candidates_f);
    best_pVecs{yr}  = candidates_pV{bi};
    best_pi{yr}     = candidates_pi{bi};
    best_fvals(yr)  = bf;
    best_source{yr} = labels{bi};

    if bf < existing_fvals(yr) - 1e-4
        nImproved = nImproved + 1;
    end
end

fprintf('Improved: %d/%d years\n', nImproved, numYears);
fprintf('Mean fval: existing=%.4f, free-init=%.4f (delta=%+.4f)\n', ...
    mean(existing_fvals), mean(best_fvals), mean(best_fvals) - mean(existing_fvals));

nExist = sum(strcmp(best_source, 'existing'));
nFwd   = sum(strcmp(best_source, 'forward'));
nBwd   = sum(strcmp(best_source, 'backward'));
fprintf('Sources: existing=%d, forward=%d, backward=%d\n\n', nExist, nFwd, nBwd);

%% ====================================================================
%  POST-PROCESS AND SAVE
%  ====================================================================
fprintf('Post-processing and saving...\n');

for yr = 1:numYears
    omeg1to5 = omega_fixed(yr, 2:end)';

    omega = zeros(ThetaTot, 1);
    omega(1) = omeg1to5(2);
    omega(2) = 1 - omeg1to5(1) - omeg1to5(2) - omeg1to5(4) - omeg1to5(5);
    omega(3) = omeg1to5(4);
    omega(4) = omeg1to5(5);
    omega(5) = omeg1to5(1);

    pVecs_out = best_pVecs{yr};
    initDists_out = best_pi{yr};

    % State ordering: ensure p33 <= p44 (mover types only)
    for theta = 1:Theta-1
        if pVecs_out(8,theta) > pVecs_out(12,theta)
            old_p = pVecs_out(:,theta);
            pVecs_out(:,theta) = [old_p(1);old_p(3);old_p(2); old_p(4);old_p(6);old_p(5); ...
                              old_p(10);old_p(12);old_p(11); old_p(7);old_p(9);old_p(8)];
            old_pi = initDists_out(:, theta);
            initDists_out(:, theta) = [old_pi(1); old_pi(2); old_pi(4); old_pi(3)];
        end
    end

    % Build TMs and ergodic
    tms = zeros(4*Theta, 4);
    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs_out(:, theta);
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
    E_mass = ergVecs(3,:) + ergVecs(4,:);
    [~, allE_idx] = max(E_mass);
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

    pVecs_out = pVecs_out(:, sortOrder);
    omega(1:Theta) = omega(sortOrder);
    initDists_out = initDists_out(:, sortOrder);

    % Rebuild after sort
    tms = zeros(4*Theta, 4);
    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs_out(:, theta);
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
    res.pVecs      = pVecs_out;
    res.omega      = omega;
    res.tms        = tms;
    res.initDists  = initDists_out;
    res.ergVecs    = ergVecs;
    res.ergVecsAll = ergVecsAll;
    res.SSR2       = best_fvals(yr);
    res.exitflag   = 1;
    res.regPenalty = regPenalty;
    res.regTerms   = regTerms;
    res.regBinding = regBinding;
    res.source     = best_source{yr};

    save(fullfile(outDir, sprintf('step3_%d.mat', yr)), 'res');
end

%% Summary table
T = table(yearsList(:), existing_fvals, best_fvals, ...
    best_fvals - existing_fvals, best_source, ...
    'VariableNames', {'Year','Sequential_fval','FreeInit_fval','Delta','Source'});
writetable(T, fullfile(outDir, 'comparison_summary.csv'));
save(fullfile(outDir, 'comparison_summary.mat'), 'T', ...
    'fwd_fvals', 'bwd_fvals', 'existing_fvals', 'best_fvals', 'best_source');

fprintf('\n=== FINAL SUMMARY ===\n');
fprintf('Improved: %d/%d\n', nImproved, numYears);
fprintf('Mean fval: sequential=%.4f -> free-init=%.4f (delta=%+.4f)\n', ...
    mean(existing_fvals), mean(best_fvals), mean(best_fvals) - mean(existing_fvals));
fprintf('Sources: existing=%d, forward=%d, backward=%d\n', nExist, nFwd, nBwd);
fprintf('Total time: %.1f minutes\n', toc(tTotal)/60);


%% ====================================================================
%  LOCAL FUNCTIONS
%  ====================================================================

function [pV_out, pi_out, fval_out] = run_year_em_fmincon( ...
        year, M_data, omeg1to5, omega, cand_pVecs, cand_pi, cfg)
% Run EM from multiple starts + fmincon refinement for one year.

Theta = cfg.Theta;
ThetaTot = cfg.ThetaTot;

nCand = numel(cand_pVecs);
em_lls   = -inf(nCand, 1);
em_pV    = cell(nCand, 1);
em_piV   = cell(nCand, 1);

for s = 1:nCand
    pV0 = cand_pVecs{s};
    pi0 = cand_pi{s};

    % Clamp to valid range
    pV0 = max(pV0, 1e-6);
    pV0 = min(pV0, 1 - 1e-6);
    pi0 = max(pi0, 1e-6);

    % All-E structural zeros
    pV0([2 5 7 8 9 11], Theta) = 0;
    pi0(3, Theta) = 0;

    % Normalize pi
    for kk = 1:Theta
        pi0(:, kk) = pi0(:, kk) / sum(pi0(:, kk));
    end

    % Ensure TM row sums < 1
    for th = 1:Theta
        for row = 1:4
            idx3 = (row-1)*3 + (1:3);
            rs_val = sum(pV0(idx3, th));
            if rs_val > 1 - 1e-6
                pV0(idx3, th) = pV0(idx3, th) * (1-1e-4) / rs_val;
            end
        end
    end
    if pV0(12, Theta) < 0.975
        pV0(12, Theta) = 0.976;
    end

    [pV_em, pi_em, ll_em, ~] = em_core_step3_freeinit( ...
        M_data, cfg.Hstates, cfg.act_idx, cfg.h4h5_groups, ...
        cfg.h1_groups, cfg.direct_pairs, ...
        Theta, ThetaTot, cfg.Nstatepath, cfg.Nactpath, ...
        cfg.All_N, omega, pV0, pi0, cfg.maxEMiter, cfg.emTol, cfg.pseudoCount);

    em_lls(s)  = ll_em(end);
    em_pV{s}   = pV_em;
    em_piV{s}  = pi_em;
end

% Pick best EM
[~, bestS] = max(em_lls);
pV_out   = em_pV{bestS};
pi_out   = em_piV{bestS};
fval_out = -em_lls(bestS);

% fmincon refinement
piParams = zeros(3*Theta, 1);
for th = 1:Theta
    piParams((th-1)*3 + (1:3)) = [pi_out(2,th); pi_out(3,th); pi_out(4,th)];
end
unks0 = [pV_out(:); piParams];
unks0 = max(cfg.lb, min(cfg.ub, unks0));

select = (1:cfg.Nactpath)';
f = @(x) ModelFitFunc_FixedOmega_AllE_FreeInit( ...
    x, select, cfg.groupNum, Theta, ThetaTot, cfg.Hstates, year, M_data, omeg1to5);

[ests, fval, ~] = optimize_from_warmstart_fi(f, unks0, cfg.lb, cfg.ub, ...
    cfg.A_ineq, cfg.B_ineq, cfg.MaxIter, cfg.MaxEvals);

if fval < fval_out
    fval_out = fval;
    pV_out = reshape(ests(1:cfg.NParTM), [12, Theta]);
    for th = 1:Theta
        piIdx = cfg.NParTM + (th-1)*3 + (1:3);
        pi234 = ests(piIdx);
        pi_out(:, th) = [1-sum(pi234); pi234(:)];
    end
    pi_out(3, Theta) = 0;
    pi_out(:, Theta) = pi_out(:, Theta) / sum(pi_out(:, Theta));
end
end


function [randPV, randPi] = make_random_start(cfg)
% Generate a random starting point for EM.
Theta = cfg.Theta;
randPV = zeros(12, Theta);
randPi = zeros(4, Theta);
for kk = 1:Theta
    if kk == Theta
        p12 = rand*0.02; p14 = rand*0.02;
        p22 = rand*0.02; p24 = rand*0.02;
        p42 = rand*0.02; p44 = 0.975 + rand*0.024;
        randPV(:,kk) = [p12;0;p14;p22;0;p24;0;0;0;p42;0;p44];
    else
        p12 = rand; p13 = (1-p12)*rand; p14 = (1-p12-p13)*rand;
        p22 = rand; p23 = (1-p22)*rand; p24 = (1-p22-p23)*rand;
        p32 = rand; p33 = (1-p32)*rand; p34 = (1-p32-p33)*rand;
        p42 = rand; p43 = (1-p42)*rand; p44 = (1-p42-p43)*rand;
        randPV(:,kk) = [p12;p13;p14;p22;p23;p24;p32;p33;p34;p42;p43;p44];
    end
    r = -log(rand(4,1));
    if kk == Theta, r(3) = 0; end
    randPi(:,kk) = r / sum(r);
end
end


function [ests, fval, ef] = optimize_from_warmstart_fi(f, x0, lb, ub, ...
    A_ineq, B_ineq, MaxIter, MaxEvals)
% 3-pass fmincon refinement for 60-param free-init objective.

opts1 = optimoptions('fmincon', ...
    'Algorithm','interior-point', 'Display','off', ...
    'MaxIterations',MaxIter, 'MaxFunctionEvaluations',MaxEvals, ...
    'OptimalityTolerance',1e-4, 'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-9, ...
    'FiniteDifferenceType','forward', 'FiniteDifferenceStepSize',1e-5);

[ests, fval, ef, out1] = fmincon(f, x0, A_ineq, B_ineq, [], [], lb, ub, [], opts1);

if ef == 2 || ef == 0 || (isfield(out1,'firstorderopt') && out1.firstorderopt > 1e-3)
    opts2 = optimoptions('fmincon', ...
        'Algorithm','sqp', 'Display','off', ...
        'MaxIterations',MaxIter, 'MaxFunctionEvaluations',MaxEvals, ...
        'OptimalityTolerance',1e-4, 'StepTolerance',1e-10, ...
        'ConstraintTolerance',1e-9, ...
        'FiniteDifferenceType','central', 'FiniteDifferenceStepSize',1e-4, ...
        'HessianApproximation','lbfgs');
    [e2, f2, ef2] = fmincon(f, ests, A_ineq, B_ineq, [], [], lb, ub, [], opts2);
    if f2 < fval
        ests = min(ub, max(lb, e2)); fval = f2; ef = ef2;
    end
end

ests_clean = min(ub, max(lb, ests));
opts3 = optimoptions('fmincon', ...
    'Algorithm','interior-point', 'Display','off', ...
    'MaxIterations',MaxIter, 'MaxFunctionEvaluations',MaxEvals, ...
    'OptimalityTolerance',1e-6, 'StepTolerance',1e-12, ...
    'ConstraintTolerance',1e-10, ...
    'FiniteDifferenceType','central', 'FiniteDifferenceStepSize',1e-6);
[e3, f3, ef3] = fmincon(f, ests_clean, A_ineq, B_ineq, [], [], lb, ub, [], opts3);
if f3 < fval + 1e-8
    ests = e3; fval = f3; ef = ef3;
end
end


function [pVecs, piVecs, ll_history, niter] = em_core_step3_freeinit( ...
    M_data, Hstates, act_idx, h4h5_groups, h1_groups, direct_pairs, ...
    Theta, ThetaTot, Nstatepath, Nactpath, ...
    All_N, omega, pVecs0, pi0, maxIter, tol, pseudoCount)
% Phase 3 EM core with free initial distribution.

pVecs  = pVecs0;
piVecs = pi0;
ll_history = nan(maxIter, 1);

for iter = 1:maxIter
    Hprobs  = zeros(Nstatepath, Theta);
    Mtildes = zeros(Nactpath, Theta);
    tms     = cell(Theta, 1);
    tm_pows = cell(Theta, 1);

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
        tms{theta} = tm;

        tp = zeros(4, 4, 10);
        tp(:,:,1) = eye(4);
        for k = 2:10
            tp(:,:,k) = tp(:,:,k-1) * tm;
        end
        tm_pows{theta} = tp;
        tm9 = tp(:,:,10);

        piVec = piVecs(:, theta);
        h = Hstates;
        P = piVec(h(:,1));

        B_vec = ones(Nstatepath, 1);
        for t = 2:4
            lin = sub2ind([4,4], h(:,t-1), h(:,t));
            B_vec = B_vec .* tm(lin);
        end

        lin_gap = sub2ind([4,4], h(:,4), h(:,5));
        C_vec = tm9(lin_gap);

        D_vec = ones(Nstatepath, 1);
        for t = 6:8
            lin = sub2ind([4,4], h(:,t-1), h(:,t));
            D_vec = D_vec .* tm(lin);
        end

        Hprobs(:, theta) = P .* B_vec .* C_vec .* D_vec;
        Mtildes(:, theta) = accumarray(act_idx, Hprobs(:, theta), [Nactpath, 1]);
    end

    MtildesNPP = [Mtildes, All_N];
    Mtilde = MtildesNPP * omega;

    valid = M_data > 0;
    ll = sum(M_data(valid) .* log(max(Mtilde(valid), realmin)));
    ll_history(iter) = ll;

    if iter > 1 && abs(ll - ll_history(iter-1)) < tol
        break;
    end

    Mtilde_safe = max(Mtilde, realmin);
    gamma_all = zeros(Nactpath, ThetaTot);
    for theta = 1:ThetaTot
        gamma_all(:, theta) = omega(theta) * MtildesNPP(:, theta) ./ Mtilde_safe;
    end

    for theta = 1:Theta
        tm = tms{theta};
        tp = tm_pows{theta};
        tm9 = tp(:,:,10);

        gap_counts = zeros(4, 4, 4, 4);
        for h4 = 1:4
            for h5 = 1:4
                denom = tm9(h4, h5);
                if denom < realmin, continue; end
                for s = 1:4
                    for sp = 1:4
                        val = 0;
                        for tau = 0:8
                            val = val + tp(h4, s, tau+1) * tp(sp, h5, 9-tau);
                        end
                        gap_counts(h4, h5, s, sp) = tm(s, sp) * val / denom;
                    end
                end
            end
        end

        w = M_data(act_idx) .* omega(theta) .* Hprobs(:, theta) ./ Mtilde_safe(act_idx);

        trans_counts = pseudoCount * ones(4, 4);
        if theta == Theta
            trans_counts(3, :) = 0;
            trans_counts(:, 3) = 0;
        end

        for dp = 1:size(direct_pairs, 1)
            t1 = direct_pairs(dp, 1);
            t2 = direct_pairs(dp, 2);
            lin = sub2ind([4,4], Hstates(:,t1), Hstates(:,t2));
            trans_counts = trans_counts + reshape(accumarray(lin, w, [16, 1]), [4, 4]);
        end

        for idx = 1:16
            [h4, h5] = ind2sub([4, 4], idx);
            mask = h4h5_groups{idx};
            w_sum = sum(w(mask));
            if w_sum > 0
                trans_counts = trans_counts + w_sum * squeeze(gap_counts(h4, h5, :, :));
            end
        end

        if theta == Theta
            trans_counts(3, :) = 0;
            trans_counts(:, 3) = 0;
        end

        for s = 1:4
            row_sum = sum(trans_counts(s, :));
            if row_sum > 0
                tm_row = trans_counts(s, :) / row_sum;
            else
                tm_row = [0 0 0 0];
            end
            pVecs((s-1)*3 + 1, theta) = tm_row(2);
            pVecs((s-1)*3 + 2, theta) = tm_row(3);
            pVecs((s-1)*3 + 3, theta) = tm_row(4);
        end

        if theta == Theta && pVecs(12, theta) < 0.975
            pVecs(12, theta) = 0.975;
            excess = pVecs(10,theta) + pVecs(12,theta) - 1;
            if excess > 0
                pVecs(10,theta) = max(0, pVecs(10,theta) - excess);
            end
        end

        pi_counts = zeros(4, 1);
        for s = 1:4
            pi_counts(s) = sum(w(h1_groups{s}));
        end
        pi_counts = pi_counts + pseudoCount;
        if theta == Theta
            pi_counts(3) = 0;
        end
        piVecs(:, theta) = pi_counts / sum(pi_counts);
    end
end

niter = iter;
ll_history = ll_history(1:niter);
end
