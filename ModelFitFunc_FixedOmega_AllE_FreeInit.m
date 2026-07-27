function fit = ModelFitFunc_FixedOmega_AllE_FreeInit(unks, select, groupNum, Theta, ThetaTot, Hstates, year, Data_full, omeg1to5)
% ModelFitFunc_FixedOmega_AllE_FreeInit — Phase 3 objective with free initial distribution
%
% Same as ModelFitFunc_FixedOmega_AllE but uses a FREE initial state
% distribution instead of the ergodic distribution.
%
% Parameters:
%   unks(1:48)  = transition matrix params (12 per type x 4 types)
%   unks(49:60) = initial distribution params (3 per type x 4 types)
%                 [pi(2), pi(3), pi(4)] for each type; pi(1) = 1 - sum
%
% All-E (type Theta) has structural zeros on state 3.
% Omegas are fixed (loaded externally).

M = Data_full(:,1);

reg_accum = 0;
Nobs = sum(M);

Nstatepath = 4^8;
Nmonth = 8;
Nact = 3;
Nactpath = Nact^Nmonth;
NParTM = 12 * Theta;

% ----- tolerant bounds guard & interior clamp -----
bTol  = 1e-10;
epsIn = 1e-12;

if any(unks < -bTol) || any(unks > 1 + bTol)
    fit = 1e12;
    return
end

x = min(1 - epsIn, max(epsIn, unks));

Mtildes = [];
pVecs = [];
ergVecAll3 = zeros(4,1);
ergAll = zeros(4, Theta);

% Fixed omegas
omega = zeros(ThetaTot,1);
omega(1,1) = omeg1to5(2,1);  % High-N
omega(2,1) = 1 - omeg1to5(1,1) - omeg1to5(2,1) - omeg1to5(4,1) - omeg1to5(5,1);  % High-U (residual)
omega(3,1) = omeg1to5(4,1);  % High-E
omega(4,1) = omeg1to5(5,1);  % All-E
omega(5,1) = omeg1to5(1,1);  % All-N

if any(~isfinite(omega)) || any(omega < 0) || abs(sum(omega)-1) > 1e-10
    fit = 1e12;
    return;
end

for theta = 1:Theta

    p12 = x(1+(theta-1)*12);  p13 = x(2+(theta-1)*12);  p14 = x(3+(theta-1)*12);
    p22 = x(4+(theta-1)*12);  p23 = x(5+(theta-1)*12);  p24 = x(6+(theta-1)*12);
    p32 = x(7+(theta-1)*12);  p33 = x(8+(theta-1)*12);  p34 = x(9+(theta-1)*12);
    p42 = x(10+(theta-1)*12); p43 = x(11+(theta-1)*12); p44 = x(12+(theta-1)*12);

    pVec = [p12; p13; p14; p22; p23; p24; p32; p33; p34; p42; p43; p44];
    pVecs = [pVecs, pVec]; %#ok<AGROW>

    tm = [1-p12-p13-p14, p12, p13, p14; ...
          1-p22-p23-p24, p22, p23, p24; ...
          1-p32-p33-p34, p32, p33, p34; ...
          1-p42-p43-p44, p42, p43, p44];

    if theta == Theta  % All-E: structural zeros on state 3
        tm = [1-p12-0-p14, p12, 0, p14; ...
              1-p22-0-p24, p22, 0, p24; ...
              0, 0, 0, 0; ...
              1-p42-0-p44, p42, 0, p44];
    end

    tms(4*theta-3:4*theta, :) = tm; %#ok<AGROW>

    % Ergodic distribution (for regularization only)
    if theta == Theta
        J = [1 2 4];
        tm3 = tm(J, J);
        pi3 = Ergodic(tm3);
        ergVec = zeros(4,1);
        ergVec(J) = pi3;
    else
        ergVec = Ergodic(tm);
    end

    ergVecAll3 = ergVecAll3 + omega(theta) * ergVec;
    ergAll(:, theta) = ergVec;

    % FREE initial distribution (KEY CHANGE from ergodic version)
    piIdx = NParTM + (theta-1)*3 + (1:3);
    pi234 = x(piIdx);
    piVec = [1 - sum(pi234); pi234(:)];
    if theta == Theta
        piVec(3) = 0;  % structural zero for All-E
        piVec(1) = 1 - piVec(2) - piVec(4);
    end

    % Validate init dist
    if any(piVec < -1e-10) || abs(sum(piVec) - 1) > 1e-8
        fit = 1e12;
        return;
    end
    piVec = max(piVec, 0);
    piVec = piVec / sum(piVec);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Activity path probabilities from the model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = zeros(Nstatepath, Nmonth);
    map = [3, 2, 1, 1];  % N, U, E, E
    for path = 1:Nstatepath
        for month = 1:Nmonth
            A(path, month) = map(Hstates(path, month));
        end
    end

    Hstates_prob = zeros(Nstatepath, 1);
    Cmid = tm^7;
    for path = 1:Nstatepath
        j1 = Hstates(path, 1);
        P = piVec(j1);  % FREE init dist (was ergVec(j1))
        B = 1;
        for t = 2:4
            B = B * tm(Hstates(path,t-1), Hstates(path,t));
        end
        C = tm(Hstates(path,4), :) * Cmid * tm(:, Hstates(path,5));
        D = 1;
        for t = 6:8
            D = D * tm(Hstates(path,t-1), Hstates(path,t));
        end
        Hstates_prob(path, 1) = P * B * C * D;
    end

    % Aggregate across hidden paths -> activity paths
    place_mults = Nact.^(7:-1:0)';
    Hact_pr = zeros(Nactpath, 1);
    for path = 1:Nstatepath
        activity_path_number = (A(path,:)-1) * place_mults + 1;
        Hact_pr(activity_path_number) = Hact_pr(activity_path_number) + Hstates_prob(path, 1);
    end
    Mtildes = [Mtildes, Hact_pr]; %#ok<AGROW>

end  % theta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularization: encourage type separation (based on ERGODIC, not init dist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m  = 0.02;    % margin
lam = 2.0;    % penalty weight
tau = 0.05;   % softplus smoothing
t2 = 1; t3 = 2; t4 = 3;

pos = @(z) tau .* (max(z./tau, 0) + log1p(exp(-abs(z./tau))));

% theta=1 (High-N) should dominate state 1 (OLF)
z = ergAll(1, t3) - ergAll(1, t2) + m;
reg_accum = reg_accum + lam * Nobs * pos(z).^2;
z = ergAll(1, t4) - ergAll(1, t2) + m;
reg_accum = reg_accum + lam * Nobs * pos(z).^2;

% theta=2 (High-U) should dominate state 2 (U)
z = ergAll(2, t2) - ergAll(2, t3) + m;
reg_accum = reg_accum + lam * Nobs * pos(z).^2;
z = ergAll(2, t4) - ergAll(2, t3) + m;
reg_accum = reg_accum + lam * Nobs * pos(z).^2;

% theta=3 (High-E) should dominate states {3,4} (E)
mass34 = @(t) ergAll(3, t) + ergAll(4, t);
z = mass34(t2) - mass34(t4) + m;
reg_accum = reg_accum + lam * Nobs * pos(z).^2;
z = mass34(t3) - mass34(t4) + m;
reg_accum = reg_accum + lam * Nobs * pos(z).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mixture likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

All_N = zeros(Nactpath, 1);
All_N(end) = 1;

ergVecAll = ergVecAll3 + [omega(ThetaTot); 0; 0; 0];
MtildesNPP = [Mtildes, All_N];

if min(omega) < 0
    fit = 1e12;
    return;
elseif min(min(MtildesNPP, [], 2), [], 1) < 0
    fit = 1e12;
    return;
end

% Stable log-likelihood via log-sum-exp
epsLL    = 1e-300;
logMt    = log(MtildesNPP + epsLL);
logOmega = log(omega(:))';
Z        = logMt + logOmega;

mx       = max(Z, [], 2);
logMix   = mx + log(sum(exp(Z - mx), 2));

fit      = -sum(M .* logMix);
fit      = fit + reg_accum;

end
