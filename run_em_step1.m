% run_em_step1.m — EM algorithm for Step 1 (free omegas), all 47 years
%
% For each year independently:
%   E-step: compute posterior type probabilities for each activity path
%   M-step: update omegas (closed form) and transition matrices (from
%           expected transition counts, handling the 9-step CPS gap)
%
% Multiple starting points per year:
%   1. fmincon best-of-both solution
%   2. Neighbor year t-1
%   3. Neighbor year t+1
%   4-5. Random
%
% After EM convergence, optionally refine with fmincon (EM+fmincon hybrid)
%
% Output: Paper2_Rewrite/Step1_FreeOmega_EM/

clc; clear;
tTotal = tic;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(fileparts(rewriteDir), ...
               'Replication Files and ResultsTablesFigures Files 2 (1)');
bestDir    = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
emDir      = fullfile(rewriteDir, 'Step1_FreeOmega_EM');
if ~exist(emDir, 'dir'), mkdir(emDir); end

addpath(origDir);

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
NParEst    = 12 * Theta;
Nest       = NParEst + ThetaTot;
yearsList  = setdiff(1976:2024, [1977 1993]);

% EM settings
maxEMiter    = 500;    % max EM iterations per start
emTol        = 1e-8;   % convergence: abs change in LL
nStarts      = 5;      % starting points per year
doFminconRef = true;   % refine EM solution with fmincon?
pseudoCount  = 1e-10;  % Dirichlet smoothing for transition counts

fprintf('=== Step 1 EM Algorithm: %d years, %d starts each ===\n', TotalYears, nStarts);

%% Build Hstates
Hstates = zeros(Nstatepath, Nmonth);
path = 1;
for j1=1:4, for j2=1:4, for j3=1:4, for j4=1:4
for j5=1:4, for j6=1:4, for j7=1:4, for j8=1:4
    Hstates(path,:) = [j1 j2 j3 j4 j5 j6 j7 j8];
    path = path+1;
end, end, end, end, end, end, end, end

%% Precompute activity path mapping (hidden path -> activity path index)
% State {1=N, 2=U, 3=Eshort, 4=Elong} -> Activity {3=N, 2=U, 1=E, 1=E}
map = [3, 2, 1, 1];
place_mults = Nact.^(7:-1:0)';
act_idx = zeros(Nstatepath, 1);
for p = 1:Nstatepath
    A = map(Hstates(p, :));
    act_idx(p) = (A - 1) * place_mults + 1;
end

%% Precompute grouping by (h4, h5) for efficient gap count accumulation
h4h5_lin = sub2ind([4,4], Hstates(:,4), Hstates(:,5));  % 65536x1, values 1:16
h4h5_groups = cell(16, 1);  % masks for each (h4,h5) pair
for idx = 1:16
    h4h5_groups{idx} = find(h4h5_lin == idx);
end

%% Direct transition observation pairs (consecutive months)
% Obs 1->2, 2->3, 3->4 (MIS 1-4), then obs 5->6, 6->7, 7->8 (MIS 13-16)
% Gap: obs 4->5 spans 9 transition steps (MIS 4 to MIS 13)
direct_pairs = [1,2; 2,3; 3,4; 5,6; 6,7; 7,8];

%% Polar type probabilities
All_E = zeros(Nactpath, 1); All_E(1) = 1;    % EEEEEEEE
All_N = zeros(Nactpath, 1); All_N(end) = 1;  % NNNNNNNN

%% Load data
fprintf('Loading data ...\n');
Data_full_TimeSeries = xlsread(fullfile(origDir, inputFileNames{groupNum}), 'J2:BD6562');

%% Load fmincon best-of-both results (starting points)
fprintf('Loading fmincon best-of-both results ...\n');
bestUnks = cell(TotalYears, 1);
bestFval = nan(TotalYears, 1);
for y = 1:TotalYears
    tmp = load(fullfile(bestDir, sprintf('step1_%d.mat', y)));
    bestUnks{y} = [tmp.res.pVecs(:); tmp.res.omega];
    bestFval(y) = tmp.res.SSR2;
end

%% Constraints for optional fmincon refinement
lb  = 1e-6 * ones(Nest, 1);
ub  = ones(Nest, 1);
Aeq = [zeros(1, NParEst), ones(1, ThetaTot)];
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

%% Pre-allocate results
EM_Fval      = nan(TotalYears, 1);
EM_FvalRef   = nan(TotalYears, 1);  % after fmincon refinement
EM_Iters     = nan(TotalYears, 1);
EM_BestStart = nan(TotalYears, 1);
EM_Seconds   = nan(TotalYears, 1);
EM_FOpt      = nan(TotalYears, 1);  % firstorderopt from fmincon refinement

%% ========== MAIN LOOP ==========
fprintf('\n--- Starting EM estimation for %d years ---\n', TotalYears);

parfor year = 1:TotalYears
    tStart = tic;
    M_data = Data_full_TimeSeries(:, year);  % 6561x1 vector of path frequencies (proportions summing to 1)
    N_total = sum(M_data);  % = 1 for frequencies
    select = (1:Nactpath)';

    % =====================================================================
    % BUILD STARTING POINTS
    % =====================================================================
    starts = cell(nStarts, 1);

    % 1. fmincon best-of-both
    starts{1} = bestUnks{year};

    % 2. Neighbor t-1
    if year > 1
        starts{2} = bestUnks{year-1};
    else
        starts{2} = bestUnks{year};  % fallback to self
    end

    % 3. Neighbor t+1
    if year < TotalYears
        starts{3} = bestUnks{year+1};
    else
        starts{3} = bestUnks{year};
    end

    % 4-5. Random
    for rs = 1:2
        rng(70000*year + rs, 'twister');
        Params = [];
        for kk = 1:Theta
            p12 = rand; p13 = (1-p12)*rand; p14 = (1-p12-p13)*rand;
            p22 = rand; p23 = (1-p22)*rand; p24 = (1-p22-p23)*rand;
            p32 = rand; p33 = (1-p32)*rand; p34 = (1-p32-p33)*rand;
            p42 = rand; p43 = (1-p42)*rand; p44 = (1-p42-p43)*rand;
            Params = [Params; p12;p13;p14;p22;p23;p24;p32;p33;p34;p42;p43;p44];
        end
        omeg = (1/ThetaTot) * ones(ThetaTot, 1);
        starts{3+rs} = max(1e-6, min(1-1e-6, [Params; omeg]));
    end

    % =====================================================================
    % RUN EM FROM EACH STARTING POINT
    % =====================================================================
    em_lls   = -inf(nStarts, 1);
    em_omegas = cell(nStarts, 1);
    em_pVecs  = cell(nStarts, 1);
    em_iters  = zeros(nStarts, 1);

    for s = 1:nStarts
        unks0 = starts{s};
        omega0 = unks0(NParEst+1:end);
        pVecs0 = reshape(unks0(1:NParEst), [12, Theta]);

        % Ensure valid starting point
        omega0 = max(omega0, 1e-6);
        omega0 = omega0 / sum(omega0);
        pVecs0 = max(pVecs0, 1e-6);
        pVecs0 = min(pVecs0, 1 - 1e-6);
        % Ensure row sums <= 1
        for theta = 1:Theta
            for row = 1:4
                idx3 = (row-1)*3 + (1:3);
                rs_val = sum(pVecs0(idx3, theta));
                if rs_val > 1 - 1e-6
                    pVecs0(idx3, theta) = pVecs0(idx3, theta) * (1-1e-4) / rs_val;
                end
            end
        end

        [omega_em, pVecs_em, ll_em, niter] = em_core( ...
            M_data, Hstates, act_idx, h4h5_groups, direct_pairs, ...
            Theta, ThetaTot, Nstatepath, Nactpath, ...
            All_E, All_N, omega0, pVecs0, maxEMiter, emTol, pseudoCount);

        em_lls(s)   = ll_em(end);
        em_omegas{s} = omega_em;
        em_pVecs{s}  = pVecs_em;
        em_iters(s)  = niter;
    end

    % Pick best EM start
    [bestLL, bestS] = max(em_lls);
    omega_best = em_omegas{bestS};
    pVecs_best = em_pVecs{bestS};

    % NLL for comparison with fmincon (fmincon minimizes NLL)
    em_nll = -bestLL;

    % =====================================================================
    % OPTIONAL: FMINCON REFINEMENT FROM EM SOLUTION
    % =====================================================================
    ref_nll = em_nll;
    ref_fopt = NaN;
    if doFminconRef
        % Pack EM solution into unks vector
        unks_em = [pVecs_best(:); omega_best];

        f = @(x) ModelFitFuncTimeSeries_NewTrMBoot_20250919_V2( ...
            x, select, groupNum, Theta, Hstates, year, M_data);

        opts_ref = optimoptions('fmincon', ...
            'Algorithm','interior-point', 'Display','off', ...
            'MaxIterations',10000, 'MaxFunctionEvaluations',20000, ...
            'OptimalityTolerance',0.001, 'StepTolerance',1e-10, ...
            'ConstraintTolerance',1e-9, 'FiniteDifferenceType','forward', ...
            'FiniteDifferenceStepSize',1e-5);

        try
            [ests_ref, fval_ref, ~, outRef] = fmincon(f, unks_em, A, B, Aeq, 1, lb, ub, [], opts_ref);
            if fval_ref < em_nll
                ref_nll = fval_ref;
                pVecs_best = reshape(ests_ref(1:NParEst), [12, Theta]);
                omega_best = ests_ref(NParEst+1:end);
            end
            if isfield(outRef, 'firstorderopt')
                ref_fopt = outRef.firstorderopt;
            end
        catch
            % fmincon failed, keep EM solution
        end
    end

    % =====================================================================
    % POST-PROCESSING (same as other scripts)
    % =====================================================================
    pVecs_out = pVecs_best;
    omega_out = omega_best;

    % State ordering: ensure p33 <= p44 (STJ vs LTJ)
    for theta = 1:Theta
        if pVecs_out(8,theta) > pVecs_out(12,theta)
            old_p = pVecs_out(:,theta);
            pVecs_out(:,theta) = [old_p(1);old_p(3);old_p(2); old_p(4);old_p(6);old_p(5); ...
                              old_p(10);old_p(12);old_p(11); old_p(7);old_p(9);old_p(8)];
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

    % Recompute ergodic after type sort
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

    % =====================================================================
    % SAVE
    % =====================================================================
    SSR2 = ref_nll;
    res = struct();
    res.results    = [];  % not needed for EM
    res.pVecs      = pVecs_out;
    res.omega      = omega_out;
    res.ergVecs    = ergVecs;
    res.ergVecsAll = ergVecsAll;
    res.SSR2       = SSR2;
    res.exitflag   = 1;
    res.em_nll     = em_nll;
    res.em_iters   = em_iters(bestS);
    res.em_bestStart = bestS;
    res.em_allLLs  = em_lls;

    parsave_em(fullfile(emDir, sprintf('step1_%d.mat', year)), res);

    % Track
    EM_Fval(year)      = SSR2;
    EM_FvalRef(year)   = ref_nll;
    EM_Iters(year)     = em_iters(bestS);
    EM_BestStart(year) = bestS;
    EM_Seconds(year)   = toc(tStart);
    EM_FOpt(year)      = ref_fopt;

    fprintf('Year %d (%d): EM_NLL=%.4f, Ref_NLL=%.4f, fmincon_NLL=%.4f, iters=%d, start=%d, %.0fs\n', ...
        year, yearsList(year), em_nll, ref_nll, bestFval(year), em_iters(bestS), bestS, toc(tStart));
end

%% ========== SUMMARY ==========
DeltaFval = bestFval - EM_Fval;  % positive = EM found better
T = table(yearsList(:), EM_Fval, bestFval, DeltaFval, EM_Iters, EM_BestStart, EM_FOpt, EM_Seconds, ...
    'VariableNames', {'Year','EM_Fval','Fmincon_Fval','Delta','EM_Iters','BestStart','FirstOrdOpt','Seconds'});
disp(T);
writetable(T, fullfile(emDir, 'convergence_em.csv'));
save(fullfile(emDir, 'convergence_em.mat'), 'T');

nBetter = sum(EM_Fval < bestFval - 1e-6);
nWorse  = sum(EM_Fval > bestFval + 1e-6);
nSame   = TotalYears - nBetter - nWorse;
fprintf('\n=== EM Step 1 Complete ===\n');
fprintf('EM better:  %d / %d years\n', nBetter, TotalYears);
fprintf('EM same:    %d / %d years\n', nSame, TotalYears);
fprintf('EM worse:   %d / %d years\n', nWorse, TotalYears);
fprintf('Total time: %.1f minutes\n', toc(tTotal)/60);


%% ========================================================================
%  EM CORE FUNCTION
%  ========================================================================
function [omega, pVecs, ll_history, niter] = em_core( ...
    M_data, Hstates, act_idx, h4h5_groups, direct_pairs, ...
    Theta, ThetaTot, Nstatepath, Nactpath, ...
    All_E, All_N, omega0, pVecs0, maxIter, tol, pseudoCount)
% EM algorithm for a single year
%
% The CPS 4-8-4 rotation means:
%   Obs 1-4: consecutive months (3 direct transitions)
%   Gap: 9 transition steps between obs 4 and obs 5
%   Obs 5-8: consecutive months (3 direct transitions)
%   Total: 6 direct + 9 gap = 15 transition steps

omega = omega0;
pVecs = pVecs0;
ll_history = nan(maxIter, 1);
N_total = sum(M_data);  % = 1 for frequencies (formulas work identically for counts or frequencies)

for iter = 1:maxIter
    % =================================================================
    % COMPUTE PATH PROBABILITIES FOR EACH MOVER TYPE
    % =================================================================
    Hprobs  = zeros(Nstatepath, Theta);
    Mtildes = zeros(Nactpath, Theta);
    tms     = cell(Theta, 1);
    ergVecsLocal = zeros(4, Theta);
    tm_pows = cell(Theta, 1);  % tm^0 through tm^9

    for theta = 1:Theta
        p = pVecs(:, theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        tms{theta} = tm;
        ergVecsLocal(:, theta) = Ergodic(tm);
        ergVec = ergVecsLocal(:, theta);

        % Precompute matrix powers for gap
        tp = zeros(4, 4, 10);  % tp(:,:,k) = tm^(k-1)
        tp(:,:,1) = eye(4);
        for k = 2:10
            tp(:,:,k) = tp(:,:,k-1) * tm;
        end
        tm_pows{theta} = tp;
        tm9 = tp(:,:,10);  % tm^9

        % --- Vectorized path probability computation ---
        h = Hstates;  % 65536 x 8

        % Initial state (ergodic)
        P = ergVec(h(:,1));

        % Direct transitions: months 1->2, 2->3, 3->4
        B_vec = ones(Nstatepath, 1);
        for t = 2:4
            lin = sub2ind([4,4], h(:,t-1), h(:,t));
            B_vec = B_vec .* tm(lin);
        end

        % Gap: tm^9(h4, h5)
        lin_gap = sub2ind([4,4], h(:,4), h(:,5));
        C_vec = tm9(lin_gap);

        % Direct transitions: months 5->6, 6->7, 7->8
        D_vec = ones(Nstatepath, 1);
        for t = 6:8
            lin = sub2ind([4,4], h(:,t-1), h(:,t));
            D_vec = D_vec .* tm(lin);
        end

        Hprobs(:, theta) = P .* B_vec .* C_vec .* D_vec;

        % Aggregate to activity paths
        Mtildes(:, theta) = accumarray(act_idx, Hprobs(:, theta), [Nactpath, 1]);
    end

    % =================================================================
    % MIXTURE PROBABILITIES AND LOG-LIKELIHOOD
    % =================================================================
    MtildesNPP = [Mtildes, All_E, All_N];
    Mtilde = MtildesNPP * omega;

    % Log-likelihood
    valid = M_data > 0;
    ll = sum(M_data(valid) .* log(max(Mtilde(valid), realmin)));
    ll_history(iter) = ll;

    % Check convergence
    if iter > 1 && abs(ll - ll_history(iter-1)) < tol
        break;
    end

    % =================================================================
    % E-STEP: Posterior type probabilities
    % =================================================================
    % gamma(k, theta) = omega(theta) * P(path_k | theta) / P(path_k)
    gamma_all = zeros(Nactpath, ThetaTot);
    Mtilde_safe = max(Mtilde, realmin);
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
    omega = omega / sum(omega);  % renormalize for safety

    % =================================================================
    % M-STEP: Update transition matrices
    % =================================================================
    for theta = 1:Theta
        tm = tms{theta};
        tp = tm_pows{theta};
        tm9 = tp(:,:,10);

        % --- Precompute gap expected transition counts ---
        % For each (h4,h5) pair, E[n(s,s') | h4,h5] over the 9 gap steps
        % = tm(s,s') / tm^9(h4,h5) * sum_{tau=0}^{8} tm^tau(h4,s) * tm^{8-tau}(s',h5)
        gap_counts = zeros(4, 4, 4, 4);  % (h4, h5, s, s')
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

        % --- Compute path weights ---
        % w(h) = M(k(h)) * omega(theta) * Hprobs(h,theta) / Mtilde(k(h))
        % (note: Mtildes cancels out when combining gamma and P(h|k,theta))
        w = M_data(act_idx) .* omega(theta) .* Hprobs(:, theta) ./ Mtilde_safe(act_idx);

        % --- Accumulate direct transition counts ---
        trans_counts = pseudoCount * ones(4, 4);  % Dirichlet smoothing
        for dp = 1:size(direct_pairs, 1)
            t1 = direct_pairs(dp, 1);
            t2 = direct_pairs(dp, 2);
            lin = sub2ind([4,4], Hstates(:,t1), Hstates(:,t2));
            trans_counts = trans_counts + reshape(accumarray(lin, w, [16, 1]), [4, 4]);
        end

        % --- Accumulate gap transition counts (grouped by h4,h5) ---
        for idx = 1:16
            [h4, h5] = ind2sub([4, 4], idx);
            mask = h4h5_groups{idx};
            w_sum = sum(w(mask));
            if w_sum > 0
                trans_counts = trans_counts + w_sum * squeeze(gap_counts(h4, h5, :, :));
            end
        end

        % --- Update transition parameters from counts ---
        for s = 1:4
            row_sum = sum(trans_counts(s, :));
            tm_row = trans_counts(s, :) / row_sum;
            pVecs((s-1)*3 + 1, theta) = tm_row(2);  % p_{s,U}
            pVecs((s-1)*3 + 2, theta) = tm_row(3);  % p_{s,Eshort}
            pVecs((s-1)*3 + 3, theta) = tm_row(4);  % p_{s,Elong}
        end
    end
end

niter = iter;
ll_history = ll_history(1:niter);
end


%% ========================================================================
%  HELPER: parallel save
%  ========================================================================
function parsave_em(outFile, res)
    save(outFile, 'res');
end
