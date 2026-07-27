function res = extract_results(ests, cfg, omega, SSR2, exitflag)
% extract_results  Parse raw estimate vector into a named struct.
%
%   res = extract_results(ests, cfg, omega, SSR2, exitflag)
%
%   Returns struct with fields:
%     tm{theta}      — 4x4 transition matrix per type
%     pVec{theta}    — 12x1 parameter vector per type
%     ergVec{theta}  — 4x1 ergodic distribution per type
%     omega          — ThetaTot x 1 type shares
%     urate(theta)   — unemployment rate by type (ergVec(2)/(1-ergVec(1)) proxy)
%     SSR2, exitflag, seed, OptTols
%     results_notReLabelled — packed column (legacy format, for comparison)
%     results_reLabelled    — packed column after type relabeling

Theta    = cfg.Theta;
ThetaTot = cfg.ThetaTot;
ppt      = cfg.paramsPerType;

% --- Extract omegas ---
if cfg.omegaEstimated
    omega = ests(end-(ThetaTot-1):end);
end

% =====================================================================
% Step 1: Extract pVecs from ests
% =====================================================================
pVecs = zeros(ppt, Theta);
for theta = 1:Theta
    base = (theta-1)*ppt;
    pVecs(:,theta) = ests(base+1:base+ppt);
end

% =====================================================================
% Step 2: Order states so that p33 <= p44 (STJ vs LTJ)
% =====================================================================
for theta = 1:Theta
    if cfg.allE_estimated && theta == Theta
        continue  % All-E has no state 3, skip
    end
    if pVecs(8,theta) > pVecs(12,theta)  % p33 > p44 → swap states 3 and 4
        old = pVecs(:,theta);
        pVecs(:,theta) = [old(1); old(3); old(2);   ...  % row1: p12, p14->p13, p13->p14
                          old(4); old(6); old(5);   ...  % row2: p22, p24->p23, p23->p24
                          old(10); old(12); old(11); ... % row3: p42->p32, p44->p33, p43->p34
                          old(7); old(9); old(8)];       % row4: p32->p42, p34->p43, p33->p44
    end
end

% =====================================================================
% Step 3: Build transition matrices and ergodic distributions
% =====================================================================
tm_cell  = cell(1, Theta);
erg_cell = cell(1, Theta);
ergMat   = zeros(4, Theta);

for theta = 1:Theta
    p = pVecs(:,theta);
    if cfg.allE_estimated && theta == Theta
        tm_cell{theta} = [1-p(1)-p(3), p(1), 0,    p(3); ...
                          1-p(4)-p(6), p(4), 0,    p(6); ...
                          0,           0,    0,    0;    ...
                          1-p(10)-p(12), p(10), 0, p(12)];
        J = [1 2 4];
        pi3 = Ergodic(tm_cell{theta}(J,J));
        ev = zeros(4,1); ev(J) = pi3;
    else
        tm_cell{theta} = [1-p(1)-p(2)-p(3),  p(1), p(2), p(3); ...
                          1-p(4)-p(5)-p(6),  p(4), p(5), p(6); ...
                          1-p(7)-p(8)-p(9),  p(7), p(8), p(9); ...
                          1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        ev = Ergodic(tm_cell{theta});
    end
    erg_cell{theta} = ev;
    ergMat(:,theta) = ev;
end

% =====================================================================
% Step 4: NOT-relabelled output (legacy packed format)
% =====================================================================
res.tm_notReLabelled  = tm_cell;
res.erg_notReLabelled = erg_cell;
res.pVecs_notReLabelled = pVecs;

results_nr = pack_results(pVecs, ergMat, omega, Theta, ThetaTot, ...
                          SSR2, exitflag, cfg);
res.results_notReLabelled = results_nr;

% =====================================================================
% Step 5: Relabel types by ergodic E mass (states 3+4), then sub-sort
% =====================================================================
E_mass = ergMat(3,:) + ergMat(4,:);
[~, I] = sort(E_mass);  % ascending: High-N, High-U, High-E, [All-E]

pVecs_rl  = pVecs(:, I);
omega_rl  = omega;
omega_rl(1:Theta) = omega(I);

% Sub-sort types 1 and 2 (High-N vs High-U) by U mass
if Theta >= 2
    U_mass_12 = [ergMat(2,I(1)), ergMat(2,I(2))];
    [~, Isub] = sort(U_mass_12);
    pVecs_rl(:,1:2) = pVecs_rl(:, Isub);
    omega_rl(1:2)   = omega_rl(Isub);
end

% Rebuild tm and erg for relabelled types
tm_rl  = cell(1, Theta);
erg_rl = cell(1, Theta);
ergMat_rl = zeros(4, Theta);
for theta = 1:Theta
    p = pVecs_rl(:,theta);
    if cfg.allE_estimated && theta == Theta
        tm_rl{theta} = [1-p(1)-p(3), p(1), 0,    p(3); ...
                        1-p(4)-p(6), p(4), 0,    p(6); ...
                        0,           0,    0,    0;    ...
                        1-p(10)-p(12), p(10), 0, p(12)];
        J = [1 2 4];
        pi3 = Ergodic(tm_rl{theta}(J,J));
        ev = zeros(4,1); ev(J) = pi3;
    else
        tm_rl{theta} = [1-p(1)-p(2)-p(3),  p(1), p(2), p(3); ...
                        1-p(4)-p(5)-p(6),  p(4), p(5), p(6); ...
                        1-p(7)-p(8)-p(9),  p(7), p(8), p(9); ...
                        1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        ev = Ergodic(tm_rl{theta});
    end
    erg_rl{theta} = ev;
    ergMat_rl(:,theta) = ev;
end

results_rl = pack_results(pVecs_rl, ergMat_rl, omega_rl, Theta, ThetaTot, ...
                          SSR2, exitflag, cfg);

% =====================================================================
% Step 6: Populate output struct
% =====================================================================
res.tm       = tm_rl;
res.ergVec   = erg_rl;
res.pVecs    = pVecs_rl;
res.omega    = omega_rl;
res.SSR2     = SSR2;
res.exitflag = exitflag;
res.seed     = cfg.seed;
res.OptTols  = cfg.OptTols;
res.results_reLabelled = results_rl;

% Unemployment rate by type: erg(state2) / (erg(state2) + erg(state3) + erg(state4))
res.urate = zeros(1, Theta);
for theta = 1:Theta
    ev = erg_rl{theta};
    denom = ev(2) + ev(3) + ev(4);
    if denom > 0
        res.urate(theta) = ev(2) / denom;
    end
end

% Aggregate ergodic distribution
ergAgg = ergMat_rl * omega_rl(1:Theta);
if cfg.allE_estimated
    ergAgg = ergAgg + [omega_rl(ThetaTot); 0; 0; 0];  % All-N adds to state 1
else
    ergAgg = ergAgg + [omega_rl(ThetaTot); 0; 0; omega_rl(ThetaTot-1)];
end
res.ergAgg = ergAgg;

end

% =========================================================================
function packed = pack_results(pVecs, ergMat, omega, Theta, ThetaTot, ...
                               SSR2, exitflag, cfg)
% Pack into the legacy column format compatible with old output files.
% Layout per column (type): 12 pVec + SSR2/0 + omega + omegaAllE/0 + omegaAllN/0
%   + OptTols/0 + seed/0 + Valid + exitflag/0 + 4 ergVec = 24 rows per type
% With appended rows for aggregate ergodic (4 rows) in the _notReLabelled version.

Valid = zeros(1, Theta);
results = [pVecs; ...
           [SSR2, zeros(1, Theta-1)]; ...
           omega(1:Theta)'; ...
           [omega(ThetaTot-1), zeros(1, Theta-1)]; ...
           [omega(ThetaTot),   zeros(1, Theta-1)]; ...
           [cfg.OptTols,       zeros(1, Theta-1)]; ...
           [cfg.seed,          zeros(1, Theta-1)]; ...
           Valid; ...
           [exitflag,          zeros(1, Theta-1)]; ...
           ergMat];

packed = results(:);
end
