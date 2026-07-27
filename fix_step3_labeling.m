% fix_step3_labeling.m — Re-sort Phase 3 results to use E-U-N criterion
%
% For years where By-E and E-U-N disagree, swaps cols 1 and 2 of pVecs
% (and corresponding omega entries) so that:
%   col 1 = High-N (max OLF), col 2 = High-U (max U), col 3 = High-E (max E)
% Then rebuilds tms, ergVecs, ergVecsAll and re-saves.

clc; clear;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
step3Dir   = fullfile(rewriteDir, 'Step3_FixedOmega');

%% Settings
numYears  = 47;
yearsList = setdiff(1976:2024, [1977 1993]);
Theta     = 4;
ThetaTot  = 5;

nFixed = 0;

for yr = 1:numYears
    fname = fullfile(step3Dir, sprintf('step3_%d.mat', yr));
    tmp = load(fname);
    res = tmp.res;
    pVecs = res.pVecs;

    % Compute ergodic distributions for mover types (cols 1-3)
    ergVecs_movers = zeros(4, 3);
    for theta = 1:3
        p = pVecs(:, theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        ergVecs_movers(:, theta) = Ergodic(tm);
    end

    E_mass = ergVecs_movers(3,:) + ergVecs_movers(4,:);
    U_mass = ergVecs_movers(2,:);

    % E-U-N criterion
    [~, highE_idx] = max(E_mass);
    remaining = setdiff(1:3, highE_idx);
    [~, rel_idx] = max(U_mass(remaining));
    highU_idx = remaining(rel_idx);
    highN_idx = setdiff(remaining, highU_idx);

    eun_order = [highN_idx, highU_idx, highE_idx];

    if isequal(eun_order, [1 2 3])
        continue;  % already correct
    end

    % Re-sort: new column order for all 4 types (col 4 = All-E stays)
    new_order = [eun_order, 4];
    pVecs_new = pVecs(:, new_order);
    omega_new = res.omega;
    omega_new(1:3) = res.omega(new_order(1:3));

    % Rebuild tms and ergVecs with new ordering
    tms_new = zeros(4*Theta, 4);
    ergVecs_new = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs_new(:, theta);
        if theta == Theta  % All-E
            tm = [1-p(1)-0-p(3), p(1), 0, p(3); ...
                  1-p(4)-0-p(6), p(4), 0, p(6); ...
                  0, 0, 0, 0; ...
                  1-p(10)-0-p(12), p(10), 0, p(12)];
            J = [1 2 4]; tm3 = tm(J,J);
            pi3 = Ergodic(tm3);
            ergVec = zeros(4,1); ergVec(J) = pi3;
        else
            tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
                  1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
                  1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
                  1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
            ergVec = Ergodic(tm);
        end
        tms_new(4*theta-3:4*theta, :) = tm;
        ergVecs_new(:, theta) = ergVec;
    end

    ergVecsAll_new = ergVecs_new * omega_new(1:Theta) + [omega_new(ThetaTot); 0; 0; 0];

    % Update and save
    res.pVecs      = pVecs_new;
    res.omega      = omega_new;
    res.tms        = tms_new;
    res.ergVecs    = ergVecs_new;
    res.ergVecsAll = ergVecsAll_new;

    save(fname, 'res');
    nFixed = nFixed + 1;
    fprintf('Fixed year %d (idx %d): swapped cols [%d,%d,%d] -> [1,2,3]\n', ...
        yearsList(yr), yr, eun_order);
end

fprintf('\nDone. Fixed %d / %d years.\n', nFixed, numYears);
