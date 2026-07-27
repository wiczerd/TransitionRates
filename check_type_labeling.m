% check_type_labeling.m — Compare type-labeling criteria across Phase 1 years
%
% For each of 47 years, loads Phase 1 results and checks:
%   (A) Do the two labeling criteria agree?
%       - "By E criterion": sort all 3 movers by E-mass, label positions 1/2/3
%         as High-N / High-U / High-E
%       - "E-U-N criterion": (1) highest E-mass → High-E, (2) among remaining
%         two, highest U → High-U, (3) remaining → High-N
%   (B) Under E-U-N criterion, does each type uniquely dominate its namesake
%       activity? i.e., is High-E the ONLY type with the highest E, High-U the
%       ONLY type with the highest U, and High-N the ONLY type with the highest N?
%       Report all violations.

clc; clear;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
phase1Dir  = fullfile(rewriteDir, 'Step1_FreeOmega_Best');

%% Settings
numYears  = 47;
yearsList = setdiff(1976:2024, [1977 1993]);

%% Load and analyze
fprintf('=== Type Labeling Comparison: Phase 1, %d years ===\n\n', numYears);

nFlips = 0;
flipYears = [];
violationYears = {};

for yr = 1:numYears
    tmp = load(fullfile(phase1Dir, sprintf('step1_%d.mat', yr)));

    % Phase 1 has 3 mover types: pVecs is 12x3
    % But the saved results are ALREADY sorted by E-U-N criterion
    % (Phase 1 code applies both sorts before saving).
    % We need the PRE-SORT ergodic distributions.
    % Since saved pVecs are already sorted, we just compute ergodics from them.
    % The saved ordering IS the E-U-N ordering: pos1=High-N, pos2=High-U, pos3=High-E.

    pVecs = tmp.res.pVecs;  % 12 x 3 (already E-U-N sorted)

    % Compute ergodic distributions for each mover type
    ergVecs = zeros(4, 3);
    for theta = 1:3
        p = pVecs(:, theta);
        tm = [1-p(1)-p(2)-p(3), p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6), p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9), p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
        ergVecs(:, theta) = Ergodic(tm);
    end

    % Activity masses (collapsing hidden states to observed activities)
    N_mass = ergVecs(1,:);             % OLF (state 1)
    U_mass = ergVecs(2,:);             % U (state 2)
    E_mass = ergVecs(3,:) + ergVecs(4,:);  % E = STJ + LTJ

    % ---- E-U-N criterion (already applied in saved data) ----
    % pos1 = High-N, pos2 = High-U, pos3 = High-E
    eun_labels = {'High-N', 'High-U', 'High-E'};
    eun_order = [1, 2, 3];  % positions

    % Verify: Step 1: pos3 should have highest E
    [~, highE_idx] = max(E_mass);
    % Step 2: among remaining, highest U → High-U
    remaining = setdiff(1:3, highE_idx);
    [~, rel_idx] = max(U_mass(remaining));
    highU_idx = remaining(rel_idx);
    % Step 3: remaining → High-N
    highN_idx = setdiff(remaining, highU_idx);

    eun_map = zeros(1,3);  % eun_map(i) = label for original position i
    eun_map(highN_idx) = 1;  % High-N
    eun_map(highU_idx) = 2;  % High-U
    eun_map(highE_idx) = 3;  % High-E

    % Since saved data is already E-U-N sorted, eun_map should be [1,2,3]
    % If not, there's an issue with the saved data

    % ---- By E criterion ----
    [~, byE_sortI] = sort(E_mass);  % ascending: pos1=lowest E, pos2, pos3=highest E
    byE_map = zeros(1,3);
    byE_map(byE_sortI(1)) = 1;  % High-N (lowest E)
    byE_map(byE_sortI(2)) = 2;  % High-U (middle E)
    byE_map(byE_sortI(3)) = 3;  % High-E (highest E)

    % ---- Compare ----
    labelsMatch = isequal(eun_map, byE_map);

    if ~labelsMatch
        nFlips = nFlips + 1;
        flipYears(end+1) = yearsList(yr); %#ok<SAGROW>

        % Find which types are swapped
        swapped = find(eun_map ~= byE_map);
        fprintf('FLIP Year %d (idx %d):\n', yearsList(yr), yr);
        for s = swapped
            eun_lbl = eun_labels{eun_map(s)};
            byE_lbl = eun_labels{byE_map(s)};
            fprintf('  Type at position %d: E-U-N calls it %s, By-E calls it %s\n', ...
                s, eun_lbl, byE_lbl);
        end
        fprintf('  E-mass: [%.4f, %.4f, %.4f]\n', E_mass);
        fprintf('  U-mass: [%.4f, %.4f, %.4f]\n', U_mass);
        fprintf('  N-mass: [%.4f, %.4f, %.4f]\n', N_mass);
        fprintf('\n');
    end

    % ---- Check dominance under E-U-N criterion ----
    % Under E-U-N: pos highN_idx = High-N, pos highU_idx = High-U, pos highE_idx = High-E
    % Check: does each type have the UNIQUE maximum in its namesake activity?

    violations = {};

    % High-E should have max E (by construction, always true)
    % But let's check anyway
    if E_mass(highE_idx) < max(E_mass(setdiff(1:3, highE_idx)))
        violations{end+1} = sprintf('High-E does NOT have max E (%.3f vs %.3f)', ...
            E_mass(highE_idx), max(E_mass(setdiff(1:3, highE_idx))));
    end

    % High-U should have max U among ALL 3 movers (not just remaining 2)
    [maxU, maxU_idx] = max(U_mass);
    if maxU_idx ~= highU_idx
        violations{end+1} = sprintf('High-U does NOT have max U among all 3 movers: High-%s has U=%.3f vs High-U U=%.3f', ...
            eun_labels{eun_map(maxU_idx)}, maxU, U_mass(highU_idx));
    end

    % High-N should have max N among ALL 3 movers
    [maxN, maxN_idx] = max(N_mass);
    if maxN_idx ~= highN_idx
        violations{end+1} = sprintf('High-N does NOT have max N among all 3 movers: High-%s has N=%.3f vs High-N N=%.3f', ...
            eun_labels{eun_map(maxN_idx)}, maxN, N_mass(highN_idx));
    end

    % Check if any type leads in MORE than one activity
    leader_E = highE_idx;  % by construction
    [~, leader_U] = max(U_mass);
    [~, leader_N] = max(N_mass);
    leaders = [leader_N, leader_U, leader_E];  % which position leads N, U, E

    for pos = 1:3
        nLeads = sum(leaders == pos);
        if nLeads > 1
            activities = {'N', 'U', 'E'};
            leading = activities(leaders == pos);
            violations{end+1} = sprintf('Position %d (%s under E-U-N) leads in %d activities: %s', ...
                pos, eun_labels{eun_map(pos)}, nLeads, strjoin(leading, ', '));
        end
        if nLeads == 0
            violations{end+1} = sprintf('Position %d (%s under E-U-N) leads in NO activity', ...
                pos, eun_labels{eun_map(pos)});
        end
    end

    if ~isempty(violations)
        violationYears{end+1} = struct('year', yearsList(yr), 'idx', yr, ...
            'violations', {violations}, ...
            'E_mass', E_mass, 'U_mass', U_mass, 'N_mass', N_mass); %#ok<SAGROW>
    end
end

%% Summary
fprintf('\n========================================\n');
fprintf('SUMMARY\n');
fprintf('========================================\n\n');

fprintf('(A) LABEL FLIPS: E-U-N vs By-E-only\n');
fprintf('    %d/%d years have different labels\n', nFlips, numYears);
if nFlips > 0
    fprintf('    Flip years: %s\n', num2str(flipYears));
end

fprintf('\n(B) DOMINANCE VIOLATIONS under E-U-N criterion\n');
fprintf('    %d/%d years have at least one violation\n', numel(violationYears), numYears);
for v = 1:numel(violationYears)
    vy = violationYears{v};
    fprintf('\n  Year %d:\n', vy.year);
    fprintf('    E-mass: [%.4f, %.4f, %.4f] (High-N, High-U, High-E)\n', vy.E_mass);
    fprintf('    U-mass: [%.4f, %.4f, %.4f]\n', vy.U_mass);
    fprintf('    N-mass: [%.4f, %.4f, %.4f]\n', vy.N_mass);
    for vi = 1:numel(vy.violations)
        fprintf('    -> %s\n', vy.violations{vi});
    end
end

fprintf('\n========================================\n');
fprintf('Done.\n');
