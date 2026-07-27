function fit = ModelFitFuncTimeSeries_FreeInitDist(unks,select,groupNum,Theta,Hstates,year,Data_full)
% Objective function for Phase 1 with FREE initial distribution.
%
% Identical to ModelFitFuncTimeSeries_NewTrMBoot_20250919_V2 except:
%   - The initial state distribution pi_theta is estimated as free parameters
%     instead of being set to the ergodic distribution of tm_theta.
%
% Parameter layout (Nest = 36 + 9 + 5 = 50):
%   unks(1:12)   = type 1 transition matrix params (p12..p44)
%   unks(13:24)  = type 2 transition matrix params
%   unks(25:36)  = type 3 transition matrix params
%   unks(37:39)  = type 1 initial dist (pi2, pi3, pi4); pi1 = 1-sum
%   unks(40:42)  = type 2 initial dist
%   unks(43:45)  = type 3 initial dist
%   unks(46:50)  = omegas (5 values, sum to 1)

M = Data_full(:,1);

MaximumLikelihood = true;
Nstatepath = 4^8;
Nmonth = 8;
Nact = 3;
Nactpath = Nact^Nmonth;
ThetaTot = Theta + 2;
NParTM = 12 * Theta;       % 36 transition matrix params
NParPi = 3 * Theta;        % 9 initial distribution params

% ----- tolerant bounds guard -----
bTol  = 1e-10;
if any(unks < -bTol) || any(unks > 1 + bTol)
    fit = 1e12;
    return
end

Mtildes = [];
omega = unks(end-(ThetaTot-1):end);

for theta = 1:Theta

    % --- Transition matrix parameters ---
    p12 = unks(1+(theta-1)*12);
    p13 = unks(2+(theta-1)*12);
    p14 = unks(3+(theta-1)*12);
    p22 = unks(4+(theta-1)*12);
    p23 = unks(5+(theta-1)*12);
    p24 = unks(6+(theta-1)*12);
    p32 = unks(7+(theta-1)*12);
    p33 = unks(8+(theta-1)*12);
    p34 = unks(9+(theta-1)*12);
    p42 = unks(10+(theta-1)*12);
    p43 = unks(11+(theta-1)*12);
    p44 = unks(12+(theta-1)*12);

    tm = [1-p12-p13-p14, p12, p13, p14; ...
          1-p22-p23-p24, p22, p23, p24; ...
          1-p32-p33-p34, p32, p33, p34; ...
          1-p42-p43-p44, p42, p43, p44];

    % --- Free initial distribution (KEY CHANGE) ---
    piIdx = NParTM + (theta-1)*3 + (1:3);
    pi2 = unks(piIdx(1));
    pi3 = unks(piIdx(2));
    pi4 = unks(piIdx(3));
    pi1 = 1 - pi2 - pi3 - pi4;
    ergVec = [pi1; pi2; pi3; pi4];

    % --- Activity mapping ---
    map = [3, 2, 1, 1];  % state -> activity: N, U, E, E
    A = zeros(Nstatepath, Nmonth);
    for path = 1:Nstatepath
        for month = 1:Nmonth
            A(path, month) = map(Hstates(path, month));
        end
    end

    % --- Path probabilities ---
    Hstates_prob = zeros(Nstatepath, 1);
    Cmid = tm^7;
    for path = 1:Nstatepath
        j1 = Hstates(path, 1);
        P = ergVec(j1);                          % FREE initial dist
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

    % --- Aggregate to activity paths ---
    place_mults = Nact.^(7:-1:0)';
    Hact_pr = zeros(Nactpath, 1);
    for path = 1:Nstatepath
        activity_path_number = (A(path,:)-1) * place_mults + 1;
        Hact_pr(activity_path_number) = Hact_pr(activity_path_number) + Hstates_prob(path);
    end
    Mtildes = [Mtildes, Hact_pr];
end

% --- Polar types ---
All_E = zeros(Nactpath, 1);
All_N = All_E;
All_E(1) = 1;
All_N(end) = 1;

MtildesNPP = [Mtildes, All_E, All_N];

% --- Log-likelihood ---
if MaximumLikelihood
    if min(omega) < 0
        disp 'Bad omega'
        fit = 1e12;
        return;
    elseif min(min(MtildesNPP,[],2),[],1) < 0
        disp 'Bad Mtilde'
        fit = 1e12;
        return;
    end

    Mtilde = MtildesNPP * omega;
    Mtilde = Mtilde(select);

    if any(~isfinite(Mtilde)) || any(Mtilde <= 0) || any(~isfinite(omega)) || any(omega < 0)
        fit = 1e12;
        return
    end

    epsLL = realmin;
    fit = -sum(M .* log(max(Mtilde, epsLL)));
else
    Mtilde = MtildesNPP * omega;
    Mtilde = Mtilde(select);
    errorM = M - Mtilde;
    fit = 0.5 * (errorM' * errorM);
end
