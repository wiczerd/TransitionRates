function candidates = build_initial_guesses(cfg, ub)
% build_initial_guesses  Generate multi-start candidate vectors for fmincon.
%
%   candidates = build_initial_guesses(cfg, ub)
%
%   Returns a cell array of column vectors, each a valid starting point.

Theta = cfg.Theta;
ppt   = cfg.paramsPerType;  % 12
s     = cfg.seed;
epsIn = 1e-6;

if cfg.omegaEstimated
    Ntotal = Theta*ppt + cfg.ThetaTot;
else
    Ntotal = Theta*ppt;
end

candidates = {};

% =========================================================================
% (a) Seeded random start — respects transition matrix structure
% =========================================================================
Params = [];
for kk = 1:Theta
    rng(s + kk, 'twister');

    if cfg.allE_estimated && kk == Theta
        % All-E: structural zeros in state 3, p44 >= 0.975
        % Row 4 constraint: p42 + p44 <= 1, p44 >= 0.975 → p42 <= 0.025
        p12 = rand*0.3; p13 = 0;       p14 = (1-p12)*rand*0.3;
        p22 = rand*0.3; p23 = 0;       p24 = (1-p22)*rand*0.3;
        p32 = 0;         p33 = 0;       p34 = 0;
        p44 = cfg.allE_p44_lb + (1-cfg.allE_p44_lb)*rand;  % [0.975, 1]
        p42 = (1-p44)*rand; p43 = 0;
    else
        p12 = rand;
        p13 = (1-p12)*rand;
        p14 = (1-p12-p13)*rand;
        p22 = rand;
        p23 = (1-p22)*rand;
        p24 = (1-p22-p23)*rand;
        p32 = rand;
        p33 = (1-p32)*rand;
        p34 = (1-p32-p33)*rand;
        p42 = rand;
        p43 = (1-p42)*rand;
        p44 = (1-p42-p43)*rand;
    end
    Params = [Params; p12;p13;p14; p22;p23;p24; p32;p33;p34; p42;p43;p44]; %#ok<AGROW>
end

unks_rand = Params;
if cfg.omegaEstimated
    unks_rand = [unks_rand; (1/cfg.ThetaTot)*ones(cfg.ThetaTot,1)];
end
unks_rand = min(1, max(0, unks_rand));
candidates{end+1} = unks_rand;

% =========================================================================
% (b) Near-uniform interior start
% =========================================================================
shrink = 0.10;
unks_u = max(epsIn, min(1-epsIn, 0.25*(1-shrink))) * ones(Ntotal, 1);

% Enforce All-E structural zeros AND row-sum feasibility
if cfg.allE_estimated
    allE_base = (Theta-1)*ppt;
    unks_u(allE_base + cfg.allE_zeroIdx) = 0;
    % p44 >= 0.975, so rows 1,2,4 have very little room for other params
    unks_u(allE_base + 12) = cfg.allE_p44_lb;       % p44 = 0.975
    % For All-E rows: each row can sum to at most 1
    % Row 1: p12 + 0 + p14 <= 1  → set both small
    % Row 2: p22 + 0 + p24 <= 1  → set both small
    % Row 4: p42 + 0 + p44 <= 1  → p42 <= 1 - p44 = 0.025
    small_val = 0.01;
    unks_u(allE_base + 1)  = small_val;  % p12
    unks_u(allE_base + 3)  = small_val;  % p14
    unks_u(allE_base + 4)  = small_val;  % p22
    unks_u(allE_base + 6)  = small_val;  % p24
    unks_u(allE_base + 10) = small_val;  % p42
end
% Enforce ub=0 constraints
unks_u = min(unks_u, ub);

if cfg.omegaEstimated
    unks_u(end-cfg.ThetaTot+1:end) = 1/cfg.ThetaTot;
end
candidates{end+1} = unks_u;

% =========================================================================
% (c) Warm-start from guess_14.mat (if available)
% =========================================================================
if isfile(cfg.guessFile)
    try
        tmp = load(cfg.guessFile);
        if isfield(tmp, 'results_notReLabelled')
            solvedGuess = reshape(tmp.results_notReLabelled, [], Theta);
            prev = solvedGuess(1:ppt, 1:Theta);
            prev_vec = prev(:);
            prev_vec = min(1, max(0, prev_vec));
            % Enforce structural zeros
            if cfg.allE_estimated
                allE_base = (Theta-1)*ppt;
                prev_vec(allE_base + cfg.allE_zeroIdx) = 0;
                prev_vec(allE_base + 12) = max(prev_vec(allE_base+12), cfg.allE_p44_lb);
            end
            prev_vec = min(prev_vec, ub(1:Theta*ppt));
            if cfg.omegaEstimated
                prev_vec = [prev_vec; (1/cfg.ThetaTot)*ones(cfg.ThetaTot,1)];
            end
            candidates{end+1} = prev_vec;
        end
    catch
        % skip if file doesn't load properly
    end
end

% =========================================================================
% (d) Jittered variants around uniform start
% =========================================================================
for j = 1:2
    jitter = (rand(Ntotal, 1) - 0.5) * 0.02;
    uJ = max(epsIn, min(1-epsIn, candidates{2} + jitter));
    if cfg.allE_estimated
        allE_base = (Theta-1)*ppt;
        uJ(allE_base + cfg.allE_zeroIdx) = 0;
        uJ(allE_base + 12) = max(uJ(allE_base+12), cfg.allE_p44_lb);
    end
    uJ = min(uJ, ub);
    if cfg.omegaEstimated
        % Renormalize omegas to sum to 1
        omIdx = Ntotal-cfg.ThetaTot+1:Ntotal;
        uJ(omIdx) = max(epsIn, uJ(omIdx));
        uJ(omIdx) = uJ(omIdx) / sum(uJ(omIdx));
    end
    candidates{end+1} = uJ; %#ok<AGROW>
end

end
