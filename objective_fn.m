function fit = objective_fn(x, cfg, Hstates, Data_full, omega)
% objective_fn  Negative log-likelihood for the hidden Markov mixture model.
%
%   fit = objective_fn(x, cfg, Hstates, Data_full, omega)
%
%   x         — parameter vector (transition matrix entries, and omegas if free mode)
%   cfg       — config struct from config_Paper2
%   Hstates   — 4^8 x 8 matrix of all hidden state paths
%   Data_full — 6561 x 1 vector of observed activity-path counts for one year
%   omega     — ThetaTot x 1 type shares (passed in for fixed mode; ignored for free mode)

Theta    = cfg.Theta;
ThetaTot = cfg.ThetaTot;
Nstates  = cfg.Nstates;
Nmonth   = cfg.Nmonth;
Nact     = cfg.Nact;
Nstatepath = cfg.Nstatepath;
Nactpath   = cfg.Nactpath;
map        = cfg.map;
ppt        = cfg.paramsPerType;  % 12

M = Data_full(:,1);

% --- Bounds guard ---
bTol  = 1e-10;
if any(x < -bTol) || any(x > 1 + bTol)
    fit = 1e12;
    return
end
epsIn = 1e-12;
x = min(1 - epsIn, max(epsIn, x));

% --- Extract omega ---
if cfg.omegaEstimated  % free_omega mode: omega is at the end of x
    omega = x(end-(ThetaTot-1):end);
end

% Guard: only check for bad omega values, NOT sum-to-1.
% fmincon enforces sum-to-1 via Aeq; its finite-difference gradient
% perturbs one omega at a time, temporarily violating the sum constraint.
if any(~isfinite(omega)) || any(omega < -epsIn)
    fit = 1e12;
    return
end

% --- Activity mapping (vectorized, computed once) ---
A = map(Hstates);  % Nstatepath x Nmonth

% --- Place multipliers for aggregating state paths -> activity paths ---
place_mults = Nact.^((Nmonth-1):-1:0)';

% --- Precompute per-type activity-path probabilities ---
Mtildes = zeros(Nactpath, Theta);
ergAll  = zeros(Nstates, Theta);
reg_accum = 0;

for theta = 1:Theta
    base = (theta-1)*ppt;

    p12 = x(base+1);  p13 = x(base+2);  p14 = x(base+3);
    p22 = x(base+4);  p23 = x(base+5);  p24 = x(base+6);
    p32 = x(base+7);  p33 = x(base+8);  p34 = x(base+9);
    p42 = x(base+10); p43 = x(base+11); p44 = x(base+12);

    % Build transition matrix
    if cfg.allE_estimated && theta == Theta
        % All-E: state 3 structurally removed
        tm = [1-p12-p14, p12, 0, p14; ...
              1-p22-p24, p22, 0, p24; ...
              0,         0,   0, 0;   ...
              1-p42-p44, p42, 0, p44];
    else
        tm = [1-p12-p13-p14, p12, p13, p14; ...
              1-p22-p23-p24, p22, p23, p24; ...
              1-p32-p33-p34, p32, p33, p34; ...
              1-p42-p43-p44, p42, p43, p44];
    end

    % Ergodic distribution
    if cfg.allE_estimated && theta == Theta
        J = [1 2 4];
        pi3 = Ergodic(tm(J,J));
        ergVec = zeros(4,1);
        ergVec(J) = pi3;
    else
        ergVec = Ergodic(tm);
    end
    ergAll(:,theta) = ergVec;

    % Path probabilities: P * B * C * D
    Cmid = tm^7;
    Hstates_prob = zeros(Nstatepath, 1);
    for path = 1:Nstatepath
        j1 = Hstates(path,1);
        P = ergVec(j1);
        B = 1;
        for t = 2:4
            B = B * tm(Hstates(path,t-1), Hstates(path,t));
        end
        C = tm(Hstates(path,4),:) * Cmid * tm(:,Hstates(path,5));
        D = 1;
        for t = 6:8
            D = D * tm(Hstates(path,t-1), Hstates(path,t));
        end
        Hstates_prob(path) = P * B * C * D;
    end

    % Aggregate to activity paths
    Hact_pr = zeros(Nactpath, 1);
    for path = 1:Nstatepath
        idx = (A(path,:)-1) * place_mults + 1;
        Hact_pr(idx) = Hact_pr(idx) + Hstates_prob(path);
    end
    Mtildes(:,theta) = Hact_pr;
end

% --- Polar types ---
if cfg.allE_estimated
    % Only All-N is hardwired (last column in omega)
    All_N = zeros(Nactpath, 1);
    All_N(end) = 1;  % NNNNNNNN is the last activity path
    MtildesAll = [Mtildes, All_N];
else
    % Both All-E and All-N are hardwired
    All_E = zeros(Nactpath, 1);
    All_N = zeros(Nactpath, 1);
    All_E(1) = 1;    % EEEEEEEE is the first activity path
    All_N(end) = 1;
    MtildesAll = [Mtildes, All_E, All_N];
end

% --- Regularization ---
if cfg.useRegularization && cfg.allE_estimated && Theta >= 3
    Nobs = sum(M);
    lam = cfg.reg.lambda;
    m   = cfg.reg.margin;
    tau = cfg.reg.tau;
    pos = @(z) tau .* (max(z./tau, 0) + log1p(exp(-abs(z./tau))));

    t2 = 1; t3 = 2; t4 = 3;  % type indices: High-N, High-U, High-E

    % Type 1 (High-N) dominates state 1
    reg_accum = reg_accum + lam*Nobs * pos(ergAll(1,t3)-ergAll(1,t2)+m).^2;
    reg_accum = reg_accum + lam*Nobs * pos(ergAll(1,t4)-ergAll(1,t2)+m).^2;

    % Type 2 (High-U) dominates state 2
    reg_accum = reg_accum + lam*Nobs * pos(ergAll(2,t2)-ergAll(2,t3)+m).^2;
    reg_accum = reg_accum + lam*Nobs * pos(ergAll(2,t4)-ergAll(2,t3)+m).^2;

    % Type 3 (High-E) dominates states {3,4}
    mass34 = @(t) ergAll(3,t) + ergAll(4,t);
    reg_accum = reg_accum + lam*Nobs * pos(mass34(t2)-mass34(t4)+m).^2;
    reg_accum = reg_accum + lam*Nobs * pos(mass34(t3)-mass34(t4)+m).^2;
end

% --- Log-likelihood (direct weighted sum, matching old code) ---
if any(omega < -epsIn)
    fit = 1e12; return
end
if any(min(MtildesAll,[],2) < -epsIn)
    fit = 1e12; return
end

Mtilde = MtildesAll * omega(:);

% Guards
if any(~isfinite(Mtilde)) || any(Mtilde <= 0) || any(~isfinite(omega)) || any(omega < 0)
    fit = 1e12; return
end

epsLL = realmin;
fit = -sum(M .* log(max(Mtilde, epsLL))) + reg_accum;

if ~isfinite(fit)
    fit = 1e12;
end
end
