% compare_year14.m — Compare new vs old output for year 14
clc; clear all;
config_Paper2;

oldFile = fullfile(cfg.oldOutputDir, 'outputFO_notReLabelled_AllE_V3_14.mat');
newFile = fullfile(cfg.outputDir, 'outputFO_AllE_14.mat');

old = load(oldFile);
new = load(newFile);

old_vec = old.results_notReLabelled;
new_vec = new.results_notReLabelled;

fprintf('Old vector length: %d\n', numel(old_vec));
fprintf('New vector length: %d\n', numel(new_vec));

% Old layout: 28 rows x 4 types = 112 (rows 1-12: pVec, 13: SSR2/0, 14-16: omegas,
% 17-18: omega polar, 19: OptTols, 20: seed, 21: valid, 22-25: ergVec, 26-28: more+aggErg)
% Actually it's 28 per type for 4 types = 112
Theta = 4;
old_mat = reshape(old_vec, [], Theta);  % should be 28 x 4
fprintf('Old reshaped: %d x %d\n', size(old_mat));

new_mat = reshape(new_vec, [], Theta);  % should be 24 x 4
fprintf('New reshaped: %d x %d\n', size(new_mat));

% Old layout: 28 rows x 4 types
% rows 1-12: pVec, 13: SSR2, 14: omega(theta), 15: omegaAllE, 16: omegaAllN,
% 17: OptTols, 18: seed, 19: Valid, 20: exitflag, 21-24: ergVec, 25-28: ergVecAll
% New layout: 24 rows x 4 types (same but no rows 25-28)

% Compare SSR2
fprintf('\nOld SSR2: %.8f\n', old_mat(13,1));
fprintf('New SSR2: %.8f\n', new_mat(13,1));
fprintf('Difference: %.6e\n', abs(old_mat(13,1) - new_mat(13,1)));

% Compare exitflags
fprintf('\nOld exitflag: %.0f\n', old_mat(20,1));
fprintf('New exitflag: %.0f\n', new_mat(20,1));

% Compare omegas
fprintf('\nOld omegas (row 14): '); fprintf('%.4f ', old_mat(14,:)); fprintf('\n');
fprintf('New omegas (row 14): '); fprintf('%.4f ', new_mat(14,:)); fprintf('\n');

% Compare pVecs (rows 1-12)
fprintf('\n--- pVec comparison (rows 1-12, not relabelled) ---\n');
for th = 1:Theta
    old_p = old_mat(1:12, th);
    new_p = new_mat(1:12, th);
    maxd = max(abs(old_p - new_p));
    fprintf('Type %d: max|diff| = %.6e\n', th, maxd);
end

% Compare ergodic distributions (rows 21-24 in both)
fprintf('\n--- Ergodic distribution comparison ---\n');
for th = 1:Theta
    old_erg = old_mat(21:24, th);
    new_erg = new_mat(21:24, th);
    maxd = max(abs(old_erg - new_erg));
    fprintf('Type %d: old=[%.4f %.4f %.4f %.4f]  new=[%.4f %.4f %.4f %.4f]  max|diff|=%.4e\n', ...
        th, old_erg, new_erg, maxd);
end

% Compare u-rates
fprintf('\n--- Unemployment rates ---\n');
for th = 1:Theta
    old_erg = old_mat(21:24, th);
    new_erg = new_mat(21:24, th);
    old_denom = old_erg(2)+old_erg(3)+old_erg(4);
    new_denom = new_erg(2)+new_erg(3)+new_erg(4);
    old_u = 0; new_u = 0;
    if old_denom > 0, old_u = old_erg(2)/old_denom; end
    if new_denom > 0, new_u = new_erg(2)/new_denom; end
    fprintf('Type %d: old_urate=%.4f  new_urate=%.4f  diff=%.4e\n', th, old_u, new_u, abs(old_u-new_u));
end
