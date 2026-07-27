% compute_fit_stats.m — Compute model fit statistics for Section 6.5
% Computes: data entropy, KL divergence, effective R², path-specific fit

addpath('.');
dataFile = 'DataTimeSeriesPM_1976_2023.xlsx';
years_all = [1976, 1978:1992, 1994:2024];
Nyears = length(years_all);

% Read data
raw = readmatrix(dataFile, 'Range', 'J2:BD6562');
Npaths = 6561;

fvals = zeros(Nyears,1);
H_data = zeros(Nyears,1);
KL_div = zeros(Nyears,1);
top5_fit = zeros(Nyears,5);  % model prob for top 5 paths

for yy = 1:Nyears
    % Load data
    M = raw(:, yy);
    M = M / sum(M);  % ensure normalized
    
    % Data entropy
    idx = M > 0;
    H_data(yy) = -sum(M(idx) .* log(M(idx)));
    
    % Load model result
    S = load(sprintf('Step3_FixedOmega_Sequential/step3_%d.mat', yy));
    fvals(yy) = S.fval;
    
    % KL divergence = cross-entropy - entropy = fval - H(data)
    KL_div(yy) = fvals(yy) - H_data(yy);
end

% Summary statistics
fprintf('\n=== Model Fit Statistics (Phase 3, Sequential) ===\n');
fprintf('Cross-entropy (fval): mean=%.4f, min=%.4f, max=%.4f\n', mean(fvals), min(fvals), max(fvals));
fprintf('Data entropy H(data): mean=%.4f, min=%.4f, max=%.4f\n', mean(H_data), min(H_data), max(H_data));
fprintf('KL divergence:        mean=%.4f, min=%.4f, max=%.4f\n', mean(KL_div), min(KL_div), max(KL_div));

% Pseudo R-squared: 1 - fval/H_uniform where H_uniform = log(6561)
H_uniform = log(Npaths);
pseudoR2 = 1 - fvals / H_uniform;
fprintf('H(uniform) = %.4f\n', H_uniform);
fprintf('Pseudo R2 (vs uniform): mean=%.4f, min=%.4f, max=%.4f\n', mean(pseudoR2), min(pseudoR2), max(pseudoR2));

% Better R2: 1 - KL/H_data (fraction of data entropy explained)
R2_kl = 1 - KL_div ./ H_data;
fprintf('R2 (1 - KL/H_data):    mean=%.4f, min=%.4f, max=%.4f\n', mean(R2_kl), min(R2_kl), max(R2_kl));

% Number of parameters
fprintf('\nParameters per year: 48\n');
fprintf('Data dimensionality: 6560 (6561 paths minus 1 for summing to 1)\n');
fprintf('Total parameters across 47 years: %d\n', 48*47);

% Year-by-year table
fprintf('\nYear  fval    H_data  KL_div  R2_kl\n');
for yy = 1:Nyears
    fprintf('%d  %.4f  %.4f  %.4f  %.4f\n', years_all(yy), fvals(yy), H_data(yy), KL_div(yy), R2_kl(yy));
end

% Save for paper
save('Step3_FixedOmega_Sequential/fit_stats.mat', 'years_all', 'fvals', 'H_data', 'KL_div', 'R2_kl', 'pseudoR2');
fprintf('\nSaved to Step3_FixedOmega_Sequential/fit_stats.mat\n');
