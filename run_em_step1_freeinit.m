% run_em_step1_freeinit.m — EM algorithm for Step 1 with FREE initial distribution
%
% Same structure as run_em_step1.m but the initial state distribution
% pi_theta is estimated as a free parameter (closed-form M-step) instead
% of being constrained to the ergodic distribution.
%
% Key advantage: the M-step for pi_theta is now EXACT (closed form),
% eliminating the main reason the previous EM underperformed.
%
% M-step for pi_theta:
%   pi_theta(s) = sum_h [w(h) * I(h1=s)] / sum_h w(h)
%
% Output: Paper2_Rewrite/Step1_FreeOmega_FreeInitDist_EM/

clc; clear;
tTotal = tic;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(fileparts(rewriteDir), ...
               'Replication Files and ResultsTablesFigures Files 2 (1)');
bestDir    = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
emDir      = fullfile(rewriteDir, 'Step1_FreeOmega_FreeInitDist_EM');
if ~exist(emDir, 'dir'), mkdir(emDir); end

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
Nact       = 3;
Nactpath   = Nact^Nmonth;
NParTM     = 12 * Theta;
NParPi     = 3 * Theta;
Nest       = NParTM + NParPi + ThetaTot;
yearsList  = setdiff(1976:2024, [1977 1993]);

% EM settings
maxEMiter    = 500;
emTol        = 1e-8;
nStarts      = 5;
doFminconRef = true;
pseudoCount  = 1e-10;

fprintf('=== Step 1 Free Init Dist (EM): %d years, %d starts each ===\n', TotalYears, nStarts);

%% Build Hstates
Hstates = zeros(Nstatepath, Nmonth);
path = 1;
for j1=1:4, for j2=1:4, for j3=1:4, for j4=1:4
for j5=1:4, for j6=1:4, for j7=1:4, for j8=1:4
    Hstates(path,:) = [j1 j2 j3 j4 j5 j6 j7 j8];
    path = path+1;
end, end, end, end, end, end, end, end

%% Precompute activity path mapping
map = [3, 2, 1, 1];
place_mults = Nact.^(7:-1:0)';
act_idx = zeros(Nstatepath, 1);
for p = 1:Nstatepath
    A = map(Hstates(p, :));
    act_idx(p) = (A - 1) * place_mults + 1;
end

%% Precompute grouping by (h4, h5)
h4h5_lin = sub2ind([4,4], Hstates(:,4), Hstates(:,5));
h4h5_groups = cell(16, 1);
for idx = 1:16
    h4h5_groups{idx} = find(h4h5_lin == idx);
end

%% Precompute grouping by h1 (for init dist M-step)
h1_groups = cell(4, 1);
for s = 1:4
    h1_groups{s} = find(Hstates(:,1) == s);
end

direct_pairs = [1,2; 2,3; 3,4; 5,6; 6,7; 7,8];

%% Polar type probabilities
All_E = zeros(Nactpath, 1); All_E(1) = 1;
All_N = zeros(Nactpath, 1); All_N(end) = 1;

%% Load data
fprintf('Loading data ...\n');
Data_full_TimeSeries = xlsread(fullfile(origDir, inputFileNames{groupNum}), 'J2:BD6562');

%% Load fmincon best-of-both results
fprintf('Loading best-of-both results ...\n');
bestUnks  = cell(TotalYears, 1);
bestFval  = nan(TotalYears, 1);
bestErgV  = cell(TotalYears, 1);
for y = 1:TotalYears
    tmp = load(fullfile(bestDir, sprintf('step1_%d.mat', y)));
    bestUnks{y} = [tmp.res.pVecs(:); tmp.res.omega];
    bestFval(y) = tmp.res.SSR2;
    ergV = zeros(4, Theta);
    for theta = 1:Theta
        p = tmp.res.pVecs(:, theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        ergV(:, theta) = Ergodic(tm);
    end
    bestErgV{y} = ergV;
end

%% Constraints for fmincon refinement
lb = 1e-6 * ones(Nest, 1);
ub = ones(Nest, 1);
Aeq = zeros(1, Nest);
Aeq(NParTM + NParPi + 1 : end) = 1;
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

%% Pre-allocate
EM_Fval      = nan(TotalYears, 1);
EM_FvalRef   = nan(TotalYears, 1);
EM_Iters     = nan(TotalYears, 1);
EM_BestStart = nan(TotalYears, 1);
EM_Seconds   = nan(TotalYears, 1);
EM_FOpt      = nan(TotalYears, 1);

%% ========== MAIN LOOP ==========
fprintf('\n--- Starting EM free-init-dist estimation for %d years ---\n', TotalYears);

parfor year = 1:TotalYears
    tStart = tic;
    M_data = Data_full_TimeSeries(:, year);
    N_total = sum(M_data);
    select = (1:Nactpath)';

    % =================================================================
    % BUILD STARTING POINTS
    % =================================================================
    starts_pVecs = cell(nStarts, 1);
    starts_omega = cell(nStarts, 1);
    starts_pi    = cell(nStarts, 1);

    % 1. fmincon best-of-both (with ergodic as starting pi)
    starts_pVecs{1} = reshape(bestUnks{year}(1:NParTM), [12, Theta]);
    starts_omega{1} = bestUnks{year}(NParTM+1:end);
    starts_pi{1}    = bestErgV{year};

    % 2. Neighbor t-1
    ny = max(year-1, 1);
    starts_pVecs{2} = reshape(bestUnks{ny}(1:NParTM), [12, Theta]);
    starts_omega{2} = bestUnks{ny}(NParTM+1:end);
    starts_pi{2}    = bestErgV{ny};

    % 3. Neighbor t+1
    ny = min(year+1, TotalYears);
    starts_pVecs{3} = reshape(bestUnks{ny}(1:NParTM), [12, Theta]);
    starts_omega{3} = bestUnks{ny}(NParTM+1:end);
    starts_pi{3}    = bestErgV{ny};

    % 4-5. Random
    for rs = 1:2
        rng(70000*year + rs, 'twister');
        Params = zeros(12, Theta);
        for kk = 1:Theta
            p12 = rand; p13 = (1-p12)*rand; p14 = (1-p12-p13)*rand;
            p22 = rand; p23 = (1-p22)*rand; p24 = (1-p22-p23)*rand;
            p32 = rand; p33 = (1-p32)*rand; p34 = (1-p32-p33)*rand;
            p42 = rand; p43 = (1-p42)*rand; p44 = (1-p42-p43)*rand;
            Params(:,kk) = [p12;p13;p14;p22;p23;p24;p32;p33;p34;p42;p43;p44];
        end
        starts_pVecs{3+rs} = Params;
        starts_omega{3+rs} = (1/ThetaTot) * ones(ThetaTot, 1);
        % Random Dirichlet(1,1,1,1) initial dist
        piRand = zeros(4, Theta);
        for kk = 1:Theta
            r = -log(rand(4,1));
            piRand(:, kk) = r / sum(r);
        end
        starts_pi{3+rs} = piRand;
    end

    % =================================================================
    % RUN EM FROM EACH STARTING POINT
    % =================================================================
    em_lls     = -inf(nStarts, 1);
    em_omegas  = cell(nStarts, 1);
    em_pVecs   = cell(nStarts, 1);
    em_piVecs  = cell(nStarts, 1);
    em_iters   = zeros(nStarts, 1);

    for s = 1:nStarts
        omega0 = starts_omega{s};
        pVecs0 = starts_pVecs{s};
        pi0    = starts_pi{s};

        % Ensure valid
        omega0 = max(omega0, 1e-6);
        omega0 = omega0 / sum(omega0);
        pVecs0 = max(pVecs0, 1e-6);
        pVecs0 = min(pVecs0, 1 - 1e-6);
        pi0 = max(pi0, 1e-6);
        for kk = 1:Theta
            pi0(:, kk) = pi0(:, kk) / sum(pi0(:, kk));
        end
        for theta = 1:Theta
            for row = 1:4
                idx3 = (row-1)*3 + (1:3);
                rs_val = sum(pVecs0(idx3, theta));
                if rs_val > 1 - 1e-6
                    pVecs0(idx3, theta) = pVecs0(idx3, theta) * (1-1e-4) / rs_val;
                end
            end
        end

        [omega_em, pVecs_em, pi_em, ll_em, niter] = em_core_freeinit( ...
            M_data, Hstates, act_idx, h4h5_groups, h1_groups, direct_pairs, ...
            Theta, ThetaTot, Nstatepath, Nactpath, ...
            All_E, All_N, omega0, pVecs0, pi0, maxEMiter, emTol, pseudoCount);

        em_lls(s)    = ll_em(end);
        em_omegas{s} = omega_em;
        em_pVecs{s}  = pVecs_em;
        em_piVecs{s} = pi_em;
        em_iters(s)  = niter;
    end

    % Pick best EM start
    [bestLL, bestS] = max(em_lls);
    omega_best = em_omegas{bestS};
    pVecs_best = em_pVecs{bestS};
    pi_best    = em_piVecs{bestS};
    em_nll = -bestLL;

    % =================================================================
    % OPTIONAL: FMINCON REFINEMENT
    % =================================================================
    ref_nll = em_nll;
    ref_fopt = NaN;
    if doFminconRef
        % Pack EM solution into unks vector [pVecs; piParams; omega]
        piParams = [];
        for theta = 1:Theta
            piParams = [piParams; pi_best(2,theta); pi_best(3,theta); pi_best(4,theta)];
        end
        unks_em = [pVecs_best(:); piParams; omega_best];

        f = @(x) ModelFitFuncTimeSeries_FreeInitDist( ...
            x, select, groupNum, Theta, Hstates, year, M_data);

        opts_ref = optimoptions('fmincon', ...
            'Algorithm','interior-point', 'Display','off', ...
            'MaxIterations',10000, 'MaxFunctionEvaluations',20000, ...
            'OptimalityTolerance',0.001, 'StepTolerance',1e-10, ...
            'ConstraintTolerance',1e-9, 'FiniteDifferenceType','forward', ...
            'FiniteDifferenceStepSize',1e-5);

        try
            [ests_ref, fval_ref, ~, outRef] = fmincon(f, unks_em, A_ineq, B_ineq, Aeq, 1, lb, ub, [], opts_ref);
            if fval_ref < em_nll
                ref_nll = fval_ref;
                pVecs_best = reshape(ests_ref(1:NParTM), [12, Theta]);
                omega_best = ests_ref(NParTM+NParPi+1:end);
                for theta = 1:Theta
                    piIdx = NParTM + (theta-1)*3 + (1:3);
                    pi234 = ests_ref(piIdx);
                    pi_best(:, theta) = [1-sum(pi234); pi234(:)];
                end
            end
            if isfield(outRef, 'firstorderopt')
                ref_fopt = outRef.firstorderopt;
            end
        catch
            % fmincon failed, keep EM solution
        end
    end

    % =================================================================
    % POST-PROCESSING
    % =================================================================
    pVecs_out = pVecs_best;
    omega_out = omega_best;
    initDists_out = pi_best;

    % State ordering: ensure p33 <= p44
    for theta = 1:Theta
        if pVecs_out(8,theta) > pVecs_out(12,theta)
            old_p = pVecs_out(:,theta);
            pVecs_out(:,theta) = [old_p(1);old_p(3);old_p(2); old_p(4);old_p(6);old_p(5); ...
                              old_p(10);old_p(12);old_p(11); old_p(7);old_p(9);old_p(8)];
            old_pi = initDists_out(:, theta);
            initDists_out(:, theta) = [old_pi(1); old_pi(2); old_pi(4); old_pi(3)];
        end
    end

    % Type ordering: sort by ergodic E mass
    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs_out(:,theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        ergVecs(:,theta) = Ergodic(tm);
    end
    E_mass = ergVecs(3,:) + ergVecs(4,:);
    [~, I] = sort(E_mass);
    pVecs_out = pVecs_out(:, I);
    omega_out(1:Theta) = omega_out(I);
    initDists_out = initDists_out(:, I);

    % Recompute ergodic
    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs_out(:,theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        ergVecs(:,theta) = Ergodic(tm);
    end

    % Sub-sort types 1 and 2 by U mass
    U_mass = ergVecs(2, 1:2);
    [~, Isub] = sort(U_mass);
    pVecs_out(:,1:2) = pVecs_out(:, Isub);
    omega_out(1:2) = omega_out(Isub);
    initDists_out(:, 1:2) = initDists_out(:, Isub);

    % Final ergodic
    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs_out(:,theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        ergVecs(:,theta) = Ergodic(tm);
    end
    ergVecsAll = ergVecs * omega_out(1:Theta) + [omega_out(ThetaTot);0;0;omega_out(ThetaTot-1)];

    % =================================================================
    % SAVE
    % =================================================================
    SSR2 = ref_nll;
    res = struct();
    res.pVecs      = pVecs_out;
    res.omega      = omega_out;
    res.initDists  = initDists_out;  % estimated initial distributions
    res.ergVecs    = ergVecs;        % ergodic (for comparison)
    res.ergVecsAll = ergVecsAll;
    res.SSR2       = SSR2;
    res.exitflag   = 1;
    res.em_nll     = em_nll;
    res.em_iters   = em_iters(bestS);
    res.em_bestStart = bestS;
    res.em_allLLs  = em_lls;

    parsave_emfi(fullfile(emDir, sprintf('step1_%d.mat', year)), res);

    % Track
    EM_Fval(year)      = SSR2;
    EM_FvalRef(year)   = ref_nll;
    EM_Iters(year)     = em_iters(bestS);
    EM_BestStart(year) = bestS;
    EM_Seconds(year)   = toc(tStart);
    EM_FOpt(year)      = ref_fopt;

    fprintf('Year %d (%d): EM_NLL=%.4f, Ref_NLL=%.4f, fmincon_best=%.4f, iters=%d, start=%d, %.0fs\n', ...
        year, yearsList(year), em_nll, ref_nll, bestFval(year), em_iters(bestS), bestS, toc(tStart));
end

%% ========== SUMMARY ==========
DeltaFval = bestFval - EM_Fval;  % positive = EM+free-init found better
T = table(yearsList(:), EM_Fval, bestFval, DeltaFval, EM_Iters, EM_BestStart, EM_FOpt, EM_Seconds, ...
    'VariableNames', {'Year','EM_FreeInit_Fval','Ergodic_Fval','Delta','EM_Iters','BestStart','FirstOrdOpt','Seconds'});
disp(T);
writetable(T, fullfile(emDir, 'convergence_em_freeinit.csv'));
save(fullfile(emDir, 'convergence_em_freeinit.mat'), 'T');

nBetter = sum(EM_Fval < bestFval - 1e-6);
nWorse  = sum(EM_Fval > bestFval + 1e-6);
nSame   = TotalYears - nBetter - nWorse;
fprintf('\n=== EM Free Init Dist Complete ===\n');
fprintf('EM+FreeInit better: %d / %d years\n', nBetter, TotalYears);
fprintf('Same:               %d / %d years\n', nSame, TotalYears);
fprintf('EM+FreeInit worse:  %d / %d years\n', nWorse, TotalYears);
fprintf('Total time: %.1f minutes\n', toc(tTotal)/60);


%% ========================================================================
%  EM CORE WITH FREE INITIAL DISTRIBUTION
%  ========================================================================
function [omega, pVecs, piVecs, ll_history, niter] = em_core_freeinit( ...
    M_data, Hstates, act_idx, h4h5_groups, h1_groups, direct_pairs, ...
    Theta, ThetaTot, Nstatepath, Nactpath, ...
    All_E, All_N, omega0, pVecs0, pi0, maxIter, tol, pseudoCount)
% EM algorithm with free initial distribution.
%
% Unlike the ergodic-constrained EM, the M-step for pi_theta is exact:
%   pi_theta(s) = sum_h [w(h) * I(h1=s)] / sum_h w(h)

omega  = omega0;
pVecs  = pVecs0;
piVecs = pi0;       % 4 x Theta matrix of initial distributions
ll_history = nan(maxIter, 1);
N_total = sum(M_data);

for iter = 1:maxIter
    % =================================================================
    % COMPUTE PATH PROBABILITIES
    % =================================================================
    Hprobs  = zeros(Nstatepath, Theta);
    Mtildes = zeros(Nactpath, Theta);
    tms     = cell(Theta, 1);
    tm_pows = cell(Theta, 1);

    for theta = 1:Theta
        p = pVecs(:, theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        tms{theta} = tm;

        % Precompute matrix powers for gap
        tp = zeros(4, 4, 10);
        tp(:,:,1) = eye(4);
        for k = 2:10
            tp(:,:,k) = tp(:,:,k-1) * tm;
        end
        tm_pows{theta} = tp;
        tm9 = tp(:,:,10);

        % FREE initial distribution (KEY CHANGE from original EM)
        piVec = piVecs(:, theta);

        h = Hstates;
        P = piVec(h(:,1));  % free init dist instead of ergodic

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

    % =================================================================
    % MIXTURE AND LOG-LIKELIHOOD
    % =================================================================
    MtildesNPP = [Mtildes, All_E, All_N];
    Mtilde = MtildesNPP * omega;

    valid = M_data > 0;
    ll = sum(M_data(valid) .* log(max(Mtilde(valid), realmin)));
    ll_history(iter) = ll;

    if iter > 1 && abs(ll - ll_history(iter-1)) < tol
        break;
    end

    % =================================================================
    % E-STEP
    % =================================================================
    Mtilde_safe = max(Mtilde, realmin);
    gamma_all = zeros(Nactpath, ThetaTot);
    for theta = 1:ThetaTot
        gamma_all(:, theta) = omega(theta) * MtildesNPP(:, theta) ./ Mtilde_safe;
    end

    % =================================================================
    % M-STEP: Update omegas
    % =================================================================
    for theta = 1:ThetaTot
        omega(theta) = sum(M_data .* gamma_all(:, theta)) / N_total;
    end
    omega = max(omega, 1e-10);
    omega = omega / sum(omega);

    % =================================================================
    % M-STEP: Update transition matrices AND initial distributions
    % =================================================================
    for theta = 1:Theta
        tm = tms{theta};
        tp = tm_pows{theta};
        tm9 = tp(:,:,10);

        % --- Gap expected transition counts ---
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

        % --- Path weights ---
        w = M_data(act_idx) .* omega(theta) .* Hprobs(:, theta) ./ Mtilde_safe(act_idx);

        % --- Direct transition counts ---
        trans_counts = pseudoCount * ones(4, 4);
        for dp = 1:size(direct_pairs, 1)
            t1 = direct_pairs(dp, 1);
            t2 = direct_pairs(dp, 2);
            lin = sub2ind([4,4], Hstates(:,t1), Hstates(:,t2));
            trans_counts = trans_counts + reshape(accumarray(lin, w, [16, 1]), [4, 4]);
        end

        % --- Gap transition counts ---
        for idx = 1:16
            [h4, h5] = ind2sub([4, 4], idx);
            mask = h4h5_groups{idx};
            w_sum = sum(w(mask));
            if w_sum > 0
                trans_counts = trans_counts + w_sum * squeeze(gap_counts(h4, h5, :, :));
            end
        end

        % --- Update transition parameters ---
        for s = 1:4
            row_sum = sum(trans_counts(s, :));
            tm_row = trans_counts(s, :) / row_sum;
            pVecs((s-1)*3 + 1, theta) = tm_row(2);
            pVecs((s-1)*3 + 2, theta) = tm_row(3);
            pVecs((s-1)*3 + 3, theta) = tm_row(4);
        end

        % --- Update initial distribution (EXACT M-step, KEY CHANGE) ---
        pi_counts = zeros(4, 1);
        for s = 1:4
            pi_counts(s) = sum(w(h1_groups{s}));
        end
        pi_counts = pi_counts + pseudoCount;  % smoothing
        piVecs(:, theta) = pi_counts / sum(pi_counts);
    end
end

niter = iter;
ll_history = ll_history(1:niter);
end


%% ========================================================================
function parsave_emfi(outFile, res)
    save(outFile, 'res');
end
