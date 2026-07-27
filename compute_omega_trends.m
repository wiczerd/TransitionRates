% compute_omega_trends.m — Phase 2: Compute trend omegas (multiple methods)
%
% Loads free omegas from Phase 1 (best-of-both, since omegas are robust)
% and fits multiple trend candidates:
%   1. Linear OLS, full sample 1976-2024
%   2. Linear OLS, 1976-2019 (pre-COVID), extrapolated
%   3. Linear OLS, excluding 2020-2021, using rest of sample
%   4. HP filter, lambda=100 (very smooth)
%   5. HP filter, lambda=6.25 (standard annual, Ravn-Uhlig)
%
% Plots: one figure per omega with actual values + all trend candidates.
% Also plots the sum of trend omegas for each method (should be ~1).

clc; clear;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
bestDir    = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
figDir     = fullfile(rewriteDir, 'Phase2_OmegaTrends', 'Figures');
outDir     = fullfile(rewriteDir, 'Phase2_OmegaTrends');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% Settings
TotalYears = 47;
ThetaTot   = 5;
yearsList  = setdiff(1976:2024, [1977 1993]);
typeNames  = {'High-N','High-U','High-E','All-E','All-N'};

%% Load omegas from best-of-both (omegas are robust across methods)
omega = nan(TotalYears, ThetaTot);
for y = 1:TotalYears
    tmp = load(fullfile(bestDir, sprintf('step1_%d.mat', y)));
    omega(y,:) = tmp.res.omega(:)';
end
fprintf('Loaded omegas for %d years from best-of-both results.\n', TotalYears);

%% Identify year indices
yr = yearsList(:);
idx_pre2020 = yr <= 2019;                       % 1976-2019 (35 years)
idx_no2021  = ~(yr == 2020 | yr == 2021);       % exclude 2020,2021 (45 years)

%% =========================================================
%  TREND 1: Linear OLS, full sample
%  =========================================================
trend1 = nan(TotalYears, ThetaTot);
for k = 1:ThetaTot
    p = polyfit(yr, omega(:,k), 1);
    trend1(:,k) = polyval(p, yr);
end

%% =========================================================
%  TREND 2: Linear OLS, 1976-2019 only, extrapolated
%  =========================================================
trend2 = nan(TotalYears, ThetaTot);
for k = 1:ThetaTot
    p = polyfit(yr(idx_pre2020), omega(idx_pre2020,k), 1);
    trend2(:,k) = polyval(p, yr);  % extrapolate to full sample
end

%% =========================================================
%  TREND 3: Linear OLS, excluding 2020-2021
%  =========================================================
trend3 = nan(TotalYears, ThetaTot);
for k = 1:ThetaTot
    p = polyfit(yr(idx_no2021), omega(idx_no2021,k), 1);
    trend3(:,k) = polyval(p, yr);
end

%% =========================================================
%  TREND 4: HP filter, lambda=100 (very smooth)
%  =========================================================
trend4 = nan(TotalYears, ThetaTot);
for k = 1:ThetaTot
    trend4(:,k) = hp_filter(omega(:,k), 100);
end

%% =========================================================
%  TREND 5: HP filter, lambda=6.25 (standard annual)
%  =========================================================
trend5 = nan(TotalYears, ThetaTot);
for k = 1:ThetaTot
    trend5(:,k) = hp_filter(omega(:,k), 6.25);
end

%% =========================================================
%  NORMALIZE: ensure trends sum to 1 for each year
%  =========================================================
trend1_n = trend1 ./ sum(trend1, 2);
trend2_n = trend2 ./ sum(trend2, 2);
trend3_n = trend3 ./ sum(trend3, 2);
trend4_n = trend4 ./ sum(trend4, 2);
trend5_n = trend5 ./ sum(trend5, 2);

%% Recession shading
recessionPeriods = [1980 1982; 1990 1991; 2001 2001; 2008 2009; 2020 2020];
addRecessionShading = @(ax) arrayfun(@(i) ...
    patch(ax, [recessionPeriods(i,1)-0.5 recessionPeriods(i,2)+0.5 ...
               recessionPeriods(i,2)+0.5 recessionPeriods(i,1)-0.5], ...
          [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5), ...
    1:size(recessionPeriods,1));

%% =========================================================
%  PLOT: One figure per omega, all trend candidates
%  =========================================================
trendNames = {'Linear (full)', 'Linear (pre-2020)', 'Linear (excl 20-21)', ...
              'HP \lambda=100', 'HP \lambda=6.25'};
trendColors = [0.8 0.2 0.2;   % red
               0.2 0.6 0.2;   % green
               0.2 0.2 0.8;   % blue
               0.8 0.5 0.0;   % orange
               0.6 0.0 0.6];  % purple

for k = 1:ThetaTot
    fig = figure('Position', [100 100 900 500]);

    % Plot actual omegas
    plot(yr, omega(:,k), 'ko-', 'MarkerSize', 4, 'LineWidth', 1.5, ...
        'MarkerFaceColor', [0.5 0.5 0.5]); hold on;

    % Plot each trend (normalized)
    plot(yr, trend1_n(:,k), '-',  'Color', trendColors(1,:), 'LineWidth', 2);
    plot(yr, trend2_n(:,k), '--', 'Color', trendColors(2,:), 'LineWidth', 2);
    plot(yr, trend3_n(:,k), ':',  'Color', trendColors(3,:), 'LineWidth', 2.5);
    plot(yr, trend4_n(:,k), '-',  'Color', trendColors(4,:), 'LineWidth', 2);
    plot(yr, trend5_n(:,k), '-',  'Color', trendColors(5,:), 'LineWidth', 2);

    title(sprintf('%s — Actual vs Trend Candidates', typeNames{k}), 'FontSize', 14);
    xlabel('Year', 'FontSize', 12);
    ylabel('Share', 'FontSize', 12);
    legend(['Actual', trendNames], 'Location', 'best', 'FontSize', 9);
    grid on;

    ax = gca;
    addRecessionShading(ax);
    ch = get(ax, 'Children'); set(ax, 'Children', flipud(ch));

    saveas(fig, fullfile(figDir, sprintf('Trend_%s.png', typeNames{k})));
end

%% =========================================================
%  PLOT: Sum of trend omegas (should be ~1 before normalization)
%  =========================================================
fig_sum = figure('Position', [100 100 900 400]);
plot(yr, sum(trend1, 2), '-',  'Color', trendColors(1,:), 'LineWidth', 2); hold on;
plot(yr, sum(trend2, 2), '--', 'Color', trendColors(2,:), 'LineWidth', 2);
plot(yr, sum(trend3, 2), ':',  'Color', trendColors(3,:), 'LineWidth', 2.5);
plot(yr, sum(trend4, 2), '-',  'Color', trendColors(4,:), 'LineWidth', 2);
plot(yr, sum(trend5, 2), '-',  'Color', trendColors(5,:), 'LineWidth', 2);
yline(1, 'k--', 'LineWidth', 1);
title('Sum of Trend Omegas (before normalization)', 'FontSize', 14);
xlabel('Year'); ylabel('Sum');
legend(trendNames, 'Location', 'best');
grid on;
saveas(fig_sum, fullfile(figDir, 'Trend_sum_check.png'));

%% =========================================================
%  PRINT SUMMARY TABLE
%  =========================================================
fprintf('\n=== TREND OMEGAS (normalized) — Selected Years ===\n');
selYears = [1976 1980 1985 1990 1995 2000 2005 2008 2010 2015 2019 2020 2024];
fprintf('\n--- Linear (full sample) ---\n');
fprintf('%6s %8s %8s %8s %8s %8s\n', 'Year', typeNames{:});
for yi = 1:length(selYears)
    idx = find(yr == selYears(yi));
    if isempty(idx), continue; end
    fprintf('%6d %8.4f %8.4f %8.4f %8.4f %8.4f\n', selYears(yi), trend1_n(idx,:));
end

fprintf('\n--- Linear (pre-2020, extrapolated) ---\n');
fprintf('%6s %8s %8s %8s %8s %8s\n', 'Year', typeNames{:});
for yi = 1:length(selYears)
    idx = find(yr == selYears(yi));
    if isempty(idx), continue; end
    fprintf('%6d %8.4f %8.4f %8.4f %8.4f %8.4f\n', selYears(yi), trend2_n(idx,:));
end

fprintf('\n--- Linear (excl 2020-2021) ---\n');
fprintf('%6s %8s %8s %8s %8s %8s\n', 'Year', typeNames{:});
for yi = 1:length(selYears)
    idx = find(yr == selYears(yi));
    if isempty(idx), continue; end
    fprintf('%6d %8.4f %8.4f %8.4f %8.4f %8.4f\n', selYears(yi), trend3_n(idx,:));
end

fprintf('\n--- HP lambda=100 ---\n');
fprintf('%6s %8s %8s %8s %8s %8s\n', 'Year', typeNames{:});
for yi = 1:length(selYears)
    idx = find(yr == selYears(yi));
    if isempty(idx), continue; end
    fprintf('%6d %8.4f %8.4f %8.4f %8.4f %8.4f\n', selYears(yi), trend4_n(idx,:));
end

%% =========================================================
%  PRINT SLOPE COMPARISON
%  =========================================================
fprintf('\n=== LINEAR SLOPES (per year, before normalization) ===\n');
fprintf('%10s %10s %10s %10s\n', '', 'Full', 'Pre-2020', 'Excl 20-21');
for k = 1:ThetaTot
    p1 = polyfit(yr, omega(:,k), 1);
    p2 = polyfit(yr(idx_pre2020), omega(idx_pre2020,k), 1);
    p3 = polyfit(yr(idx_no2021), omega(idx_no2021,k), 1);
    fprintf('%10s %+10.6f %+10.6f %+10.6f\n', typeNames{k}, p1(1), p2(1), p3(1));
end

%% Save all trends
save(fullfile(outDir, 'omega_trends.mat'), ...
    'yearsList', 'omega', 'trend1_n', 'trend2_n', 'trend3_n', 'trend4_n', 'trend5_n', ...
    'typeNames', 'trendNames');
fprintf('\nTrends saved to %s\n', fullfile(outDir, 'omega_trends.mat'));
fprintf('Figures saved to %s\n', figDir);


%% =========================================================
%  HP FILTER (local function)
%  =========================================================
function tau = hp_filter(y, lambda)
% One-sided HP filter via matrix formulation.
% tau = argmin sum((y_t - tau_t)^2) + lambda * sum((tau_{t+1} - 2*tau_t + tau_{t-1})^2)

    T = length(y);
    % Build second-difference matrix
    e = ones(T, 1);
    D = spdiags([e -2*e e], 0:2, T-2, T);
    % Solve (I + lambda * D'D) * tau = y
    tau = (speye(T) + lambda * (D' * D)) \ y(:);
    tau = tau(:);
end
