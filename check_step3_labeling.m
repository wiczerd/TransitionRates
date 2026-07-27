% check_step3_labeling.m — Diagnose Phase 3 type labeling
%
% Phase 3 results were saved with By-E post-sort (ascending E-mass).
% This script:
%   (A) Checks whether By-E ordering matches E-U-N criterion
%   (B) Under E-U-N criterion, checks dominance violations:
%       does each type uniquely dominate its namesake activity?
%   (C) Reports all findings — does NOT modify any files
%
% This mirrors check_type_labeling.m (Phase 1) but applied to Phase 3.

clc; clear;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
step3Dir   = fullfile(rewriteDir, 'Step3_FixedOmega');

%% Settings
numYears  = 47;
yearsList = setdiff(1976:2024, [1977 1993]);
Theta     = 4;  % High-N, High-U, High-E, All-E

%% Analyze each year
fprintf('=== Phase 3 Type Labeling Diagnosis ===\n\n');

nFlips = 0;
flipYears = [];
violationYears = {};

for yr = 1:numYears
    fname = fullfile(step3Dir, sprintf('step3_%d.mat', yr));
    tmp = load(fname);
    res = tmp.res;

    pVecs = res.pVecs;  % 12 x 4 (currently By-E sorted)

    % Compute ergodic distributions for all 4 types
    ergVecs = zeros(4, Theta);
    for theta = 1:Theta
        p = pVecs(:, theta);
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
        ergVecs(:, theta) = ergVec;
    end

    % Activity masses for mover types only (columns 1-3)
    N_mass = ergVecs(1, 1:3);                      % OLF
    U_mass = ergVecs(2, 1:3);                      % U
    E_mass = ergVecs(3, 1:3) + ergVecs(4, 1:3);   % E = STJ + LTJ

    % ---- (A) E-U-N criterion ----
    % Step 1: highest E among movers -> High-E
    [~, highE_idx] = max(E_mass);
    % Step 2: among remaining two, highest U -> High-U
    remaining = setdiff(1:3, highE_idx);
    [~, rel_idx] = max(U_mass(remaining));
    highU_idx = remaining(rel_idx);
    % Step 3: remaining -> High-N
    highN_idx = setdiff(remaining, highU_idx);

    eun_order = [highN_idx, highU_idx, highE_idx];
    eun_labels = {'High-N', 'High-U', 'High-E'};

    % Build eun_map: eun_map(i) = label index for current column i
    eun_map = zeros(1,3);
    eun_map(highN_idx) = 1;  % High-N
    eun_map(highU_idx) = 2;  % High-U
    eun_map(highE_idx) = 3;  % High-E

    % By-E says [1,2,3] = [High-N, High-U, High-E]
    byE_map = [1, 2, 3];

    labelsMatch = isequal(eun_map, byE_map);

    if ~labelsMatch
        nFlips = nFlips + 1;
        flipYears(end+1) = yearsList(yr); %#ok<SAGROW>

        swapped = find(eun_map ~= byE_map);
        fprintf('FLIP Year %d (idx %d):\n', yearsList(yr), yr);
        for s = swapped
            fprintf('  Col %d: By-E calls it %s, E-U-N calls it %s\n', ...
                s, eun_labels{byE_map(s)}, eun_labels{eun_map(s)});
        end
        fprintf('  E-mass: [%.4f, %.4f, %.4f]\n', E_mass);
        fprintf('  U-mass: [%.4f, %.4f, %.4f]\n', U_mass);
        fprintf('  N-mass: [%.4f, %.4f, %.4f]\n', N_mass);
        fprintf('\n');
    end

    % ---- (B) Dominance violations under E-U-N ----
    violations = {};

    % High-E should have max E (by construction of E-U-N, always true)
    if E_mass(highE_idx) < max(E_mass(setdiff(1:3, highE_idx)))
        violations{end+1} = sprintf('High-E does NOT have max E (%.4f vs %.4f)', ...
            E_mass(highE_idx), max(E_mass(setdiff(1:3, highE_idx))));
    end

    % High-U should have max U among ALL 3 movers (not guaranteed by E-U-N)
    [maxU, maxU_idx] = max(U_mass);
    if maxU_idx ~= highU_idx
        violations{end+1} = sprintf('High-U does NOT have max U among all 3: col %d (%s) has U=%.4f vs High-U (col %d) U=%.4f', ...
            maxU_idx, eun_labels{eun_map(maxU_idx)}, maxU, highU_idx, U_mass(highU_idx));
    end

    % High-N should have max N among ALL 3 movers (not guaranteed by E-U-N)
    [maxN, maxN_idx] = max(N_mass);
    if maxN_idx ~= highN_idx
        violations{end+1} = sprintf('High-N does NOT have max N among all 3: col %d (%s) has N=%.4f vs High-N (col %d) N=%.4f', ...
            maxN_idx, eun_labels{eun_map(maxN_idx)}, maxN, highN_idx, N_mass(highN_idx));
    end

    % Check if any type leads in more than one activity or zero activities
    [~, leader_E] = max(E_mass);
    [~, leader_U] = max(U_mass);
    [~, leader_N] = max(N_mass);
    leaders = [leader_N, leader_U, leader_E];  % which col leads N, U, E

    for pos = 1:3
        nLeads = sum(leaders == pos);
        if nLeads > 1
            activities = {'N', 'U', 'E'};
            leading = activities(leaders == pos);
            violations{end+1} = sprintf('Col %d (%s under E-U-N) leads in %d activities: %s', ...
                pos, eun_labels{eun_map(pos)}, nLeads, strjoin(leading, ', '));
        end
        if nLeads == 0
            violations{end+1} = sprintf('Col %d (%s under E-U-N) leads in NO activity', ...
                pos, eun_labels{eun_map(pos)});
        end
    end

    if ~isempty(violations)
        violationYears{end+1} = struct('year', yearsList(yr), 'idx', yr, ...
            'violations', {violations}, ...
            'E_mass', E_mass, 'U_mass', U_mass, 'N_mass', N_mass, ...
            'eun_order', eun_order); %#ok<SAGROW>
    end
end

%% Summary
fprintf('\n========================================\n');
fprintf('SUMMARY\n');
fprintf('========================================\n\n');

fprintf('(A) LABEL FLIPS: E-U-N vs By-E\n');
fprintf('    %d / %d years have different labels\n', nFlips, numYears);
if nFlips > 0
    fprintf('    Flip years: %s\n', num2str(flipYears));
end

fprintf('\n(B) DOMINANCE VIOLATIONS under E-U-N criterion\n');
fprintf('    %d / %d years have at least one violation\n', numel(violationYears), numYears);
for v = 1:numel(violationYears)
    vy = violationYears{v};
    fprintf('\n  Year %d (idx %d):\n', vy.year, vy.idx);
    fprintf('    E-U-N order: [%d, %d, %d] = [High-N, High-U, High-E]\n', vy.eun_order);
    fprintf('    E-mass: [%.4f, %.4f, %.4f] (current cols 1,2,3)\n', vy.E_mass);
    fprintf('    U-mass: [%.4f, %.4f, %.4f]\n', vy.U_mass);
    fprintf('    N-mass: [%.4f, %.4f, %.4f]\n', vy.N_mass);
    for vi = 1:numel(vy.violations)
        fprintf('    -> %s\n', vy.violations{vi});
    end
end

fprintf('\n========================================\n');
fprintf('Done. This script is DIAGNOSTIC ONLY — no files were modified.\n');
