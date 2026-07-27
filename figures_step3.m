% figures_step3.m — Generate Phase 3 figures for the paper
%
% Figures:
%   1. Type shares: Phase 1 free estimates + Phase 2 trends
%   2. Ergodic E-mass (employment) by mover type over time
%   3. Ergodic U-mass (unemployment) by mover type over time
%   4. Ergodic N-mass (OLF) by mover type over time
%   5. Unemployment rate [U/(U+E)] by type over time — the key "brunt" figure
%   6. Aggregate unemployment rate: model-implied vs simple aggregate

clc; clear; close all;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
step1Dir   = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
step3Dir   = fullfile(rewriteDir, 'Step3_FixedOmega');
figDir     = fullfile(step3Dir, 'Figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% Settings
numYears  = 47;
yearsList = setdiff(1976:2024, [1977 1993]);

% NBER recession years (approximate for annual data)
recessions = {[1980 1980], [1981 1982], [1990 1991], [2001 2001], ...
              [2008 2009], [2020 2020]};

% Colors
cHighN = [0.2 0.6 0.2];   % green
cHighU = [0.8 0.2 0.2];   % red
cHighE = [0.2 0.2 0.8];   % blue
cAllE  = [0.5 0.5 0.5];   % gray
cAllN  = [0.7 0.5 0.2];   % brown

%% Load omega_fixed (Phase 2 trends)
omf = load(fullfile(rewriteDir, 'omega_fixed.mat'));
omega_fixed = omf.omega_fixed;  % 47 x 6: [year, All-N, High-N, High-U, High-E, All-E]

%% Load Phase 1 and Phase 3 results
% Phase 1: omega is [High-N, High-U, High-E, All-E, All-N]
% Phase 3: omega is [High-N, High-U, High-E, All-E, All-N] (after E-U-N fix)
%          ergVecs is 4x4: [OLF; U; STJ; LTJ] x [High-N, High-U, High-E, All-E]

p1_omega = zeros(numYears, 5);   % [High-N, High-U, High-E, All-E, All-N]
p3_ergVecs = zeros(4, 4, numYears);  % 4 states x 4 types x 47 years
p3_omega   = zeros(numYears, 5);

for yr = 1:numYears
    % Phase 1
    tmp1 = load(fullfile(step1Dir, sprintf('step1_%d.mat', yr)));
    p1_omega(yr, :) = tmp1.res.omega';

    % Phase 3
    tmp3 = load(fullfile(step3Dir, sprintf('step3_%d.mat', yr)));
    p3_ergVecs(:, :, yr) = tmp3.res.ergVecs;
    p3_omega(yr, :) = tmp3.res.omega';
end

%% Compute activity masses from Phase 3 ergodic distributions
% For mover types (1-3) and All-E (4)
E_mass = squeeze(p3_ergVecs(3,:,:) + p3_ergVecs(4,:,:));  % 4 x 47
U_mass = squeeze(p3_ergVecs(2,:,:));                       % 4 x 47
N_mass = squeeze(p3_ergVecs(1,:,:));                       % 4 x 47

% Unemployment rate: U/(U+E) for each type
Urate = U_mass ./ (U_mass + E_mass);  % 4 x 47

% Aggregate unemployment rate (model-implied)
% agg_U = sum over types of omega * U_mass_type
agg_U_mass = zeros(1, numYears);
agg_E_mass = zeros(1, numYears);
for yr = 1:numYears
    om = p3_omega(yr, :);  % [High-N, High-U, High-E, All-E, All-N]
    % U mass contribution: mover types + All-E + All-N(=0)
    agg_U_mass(yr) = om(1)*U_mass(1,yr) + om(2)*U_mass(2,yr) + ...
                     om(3)*U_mass(3,yr) + om(4)*U_mass(4,yr);
    agg_E_mass(yr) = om(1)*E_mass(1,yr) + om(2)*E_mass(2,yr) + ...
                     om(3)*E_mass(3,yr) + om(4)*E_mass(4,yr);
end
agg_Urate = agg_U_mass ./ (agg_U_mass + agg_E_mass);

%% Helper function for recession shading
add_recessions = @(ax, recs, yrs) cellfun(@(r) ...
    patch(ax, [r(1)-0.5, r(2)+0.5, r(2)+0.5, r(1)-0.5], ...
    [ax.YLim(1), ax.YLim(1), ax.YLim(2), ax.YLim(2)], ...
    [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5), recs);

% =====================================================================
%% FIGURE 1: Type Shares — Phase 1 estimates + Phase 2 trends
% =====================================================================
fig1 = figure('Position', [50 50 900 600], 'Color', 'w');

typeNames = {'High-N', 'High-U', 'High-E', 'All-E', 'All-N'};
% Phase 1 omega order: [High-N, High-U, High-E, All-E, All-N]
% omega_fixed order:   [year, All-N, High-N, High-U, High-E, All-E]
% Map omega_fixed columns to Phase 1 order: [3, 4, 5, 6, 2]
trendCols = [3, 4, 5, 6, 2];
typeColors = {cHighN, cHighU, cHighE, cAllE, cAllN};

for i = 1:5
    subplot(2, 3, i);
    hold on;

    % Recession shading (set ylim first)
    yl = [0, max(p1_omega(:,i))*1.3 + 0.01];
    ylim(yl);
    for r = 1:numel(recessions)
        rec = recessions{r};
        patch([rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
              [yl(1), yl(1), yl(2), yl(2)], ...
              [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5);
    end

    % Phase 1 free estimates
    plot(yearsList, p1_omega(:,i), 'o', 'Color', typeColors{i}, ...
         'MarkerSize', 4, 'MarkerFaceColor', typeColors{i});
    % Phase 2 trend
    plot(yearsList, omega_fixed(:, trendCols(i)), '-', 'Color', 'k', ...
         'LineWidth', 2);

    title(typeNames{i}, 'FontSize', 12);
    xlim([1975 2025]);
    if i <= 3
        ylabel('Population share');
    end
    grid on; box on;
end

sgtitle('Figure 1: Type Shares — Phase 1 Estimates and Phase 2 Trends', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig1, fullfile(figDir, 'fig1_omega_trends.png'));
saveas(fig1, fullfile(figDir, 'fig1_omega_trends.fig'));

% =====================================================================
%% FIGURE 2: Ergodic Employment Mass by Mover Type
% =====================================================================
fig2 = figure('Position', [50 50 800 450], 'Color', 'w');
hold on;

yl = [0 1.05];
ylim(yl);
for r = 1:numel(recessions)
    rec = recessions{r};
    patch([rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
          [yl(1), yl(1), yl(2), yl(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5);
end

plot(yearsList, E_mass(1,:), '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, E_mass(2,:), '-o', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, E_mass(3,:), '-o', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
plot(yearsList, E_mass(4,:), '-o', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');

legend('Location', 'southwest', 'FontSize', 10);
xlabel('Year'); ylabel('Ergodic E-mass (STJ + LTJ)');
title('Figure 2: Employment Mass by Type (Phase 3)', 'FontSize', 14);
xlim([1975 2025]); grid on; box on;
saveas(fig2, fullfile(figDir, 'fig2_emass_bytype.png'));
saveas(fig2, fullfile(figDir, 'fig2_emass_bytype.fig'));

% =====================================================================
%% FIGURE 3: Ergodic Unemployment Mass by Mover Type
% =====================================================================
fig3 = figure('Position', [50 50 800 450], 'Color', 'w');
hold on;

yl = [0 max(U_mass(:))*1.15];
ylim(yl);
for r = 1:numel(recessions)
    rec = recessions{r};
    patch([rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
          [yl(1), yl(1), yl(2), yl(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5);
end

plot(yearsList, U_mass(1,:), '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, U_mass(2,:), '-o', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, U_mass(3,:), '-o', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
plot(yearsList, U_mass(4,:), '-o', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');

legend('Location', 'northeast', 'FontSize', 10);
xlabel('Year'); ylabel('Ergodic U-mass');
title('Figure 3: Unemployment Mass by Type (Phase 3)', 'FontSize', 14);
xlim([1975 2025]); grid on; box on;
saveas(fig3, fullfile(figDir, 'fig3_umass_bytype.png'));
saveas(fig3, fullfile(figDir, 'fig3_umass_bytype.fig'));

% =====================================================================
%% FIGURE 4: Ergodic OLF Mass by Mover Type
% =====================================================================
fig4 = figure('Position', [50 50 800 450], 'Color', 'w');
hold on;

yl = [0 max(N_mass(:))*1.15];
ylim(yl);
for r = 1:numel(recessions)
    rec = recessions{r};
    patch([rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
          [yl(1), yl(1), yl(2), yl(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5);
end

plot(yearsList, N_mass(1,:), '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, N_mass(2,:), '-o', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, N_mass(3,:), '-o', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
plot(yearsList, N_mass(4,:), '-o', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');

legend('Location', 'northeast', 'FontSize', 10);
xlabel('Year'); ylabel('Ergodic OLF-mass');
title('Figure 4: OLF Mass by Type (Phase 3)', 'FontSize', 14);
xlim([1975 2025]); grid on; box on;
saveas(fig4, fullfile(figDir, 'fig4_nmass_bytype.png'));
saveas(fig4, fullfile(figDir, 'fig4_nmass_bytype.fig'));

% =====================================================================
%% FIGURE 5: Unemployment Rate [U/(U+E)] by Type — THE KEY FIGURE
% =====================================================================
fig5 = figure('Position', [50 50 900 500], 'Color', 'w');
hold on;

yl = [0 min(1, max(Urate(:))*1.1)];
ylim(yl);
for r = 1:numel(recessions)
    rec = recessions{r};
    patch([rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
          [yl(1), yl(1), yl(2), yl(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5);
end

plot(yearsList, Urate(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 2, ...
     'MarkerSize', 4, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, Urate(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 2, ...
     'MarkerSize', 4, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, Urate(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 2, ...
     'MarkerSize', 4, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
plot(yearsList, Urate(4,:)*100, '-v', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');
plot(yearsList, agg_Urate*100, '-k', 'LineWidth', 2.5, 'DisplayName', 'Model aggregate');

% Fix ylim for percentage
ylim([0 max(Urate(:)*100)*1.1]);

legend('Location', 'northeast', 'FontSize', 11);
xlabel('Year', 'FontSize', 12);
ylabel('Unemployment Rate, %', 'FontSize', 12);
title('Figure 5: Unemployment Rate by Type (Phase 3)', 'FontSize', 14);
xlim([1975 2025]); grid on; box on;
set(gca, 'FontSize', 11);
saveas(fig5, fullfile(figDir, 'fig5_urate_bytype.png'));
saveas(fig5, fullfile(figDir, 'fig5_urate_bytype.fig'));

% =====================================================================
%% FIGURE 6: Unemployment Rate — Mover Types Only (zoomed in)
% =====================================================================
fig6 = figure('Position', [50 50 900 500], 'Color', 'w');
hold on;

yl = [0 max(Urate(2,:)*100)*1.1];
ylim(yl);
for r = 1:numel(recessions)
    rec = recessions{r};
    patch([rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
          [yl(1), yl(1), yl(2), yl(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5);
end

plot(yearsList, Urate(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 2, ...
     'MarkerSize', 4, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, Urate(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 2, ...
     'MarkerSize', 4, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, Urate(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 2, ...
     'MarkerSize', 4, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');

legend('Location', 'northeast', 'FontSize', 11);
xlabel('Year', 'FontSize', 12);
ylabel('Unemployment Rate, %', 'FontSize', 12);
title('Figure 6: Unemployment Rate — Mover Types (Phase 3)', 'FontSize', 14);
xlim([1975 2025]); grid on; box on;
set(gca, 'FontSize', 11);
saveas(fig6, fullfile(figDir, 'fig6_urate_movers.png'));
saveas(fig6, fullfile(figDir, 'fig6_urate_movers.fig'));

% =====================================================================
%% FIGURE 7: Ergodic Distributions — Stacked bars for selected years
% =====================================================================
selectedYears = [1978, 1982, 1997, 2007, 2009, 2019, 2020, 2024];
selIdx = arrayfun(@(y) find(yearsList == y), selectedYears);

fig7 = figure('Position', [50 50 1000 400], 'Color', 'w');
typeLabels = {'High-N', 'High-U', 'High-E', 'All-E'};
actLabels  = {'OLF', 'U', 'E'};

for ti = 1:4
    subplot(1, 4, ti);
    % For each selected year, stack [N, U, E]
    barData = zeros(numel(selectedYears), 3);
    for si = 1:numel(selectedYears)
        barData(si, :) = [N_mass(ti, selIdx(si)), ...
                          U_mass(ti, selIdx(si)), ...
                          E_mass(ti, selIdx(si))];
    end
    b = bar(barData, 'stacked');
    b(1).FaceColor = cHighN;  % OLF = green
    b(2).FaceColor = cHighU;  % U = red
    b(3).FaceColor = cHighE;  % E = blue

    set(gca, 'XTickLabel', selectedYears, 'XTickLabelRotation', 45);
    ylim([0 1.05]);
    title(typeLabels{ti}, 'FontSize', 12);
    if ti == 1
        ylabel('Ergodic share');
        legend(actLabels, 'Location', 'northwest', 'FontSize', 8);
    end
    grid on; box on;
end

sgtitle('Figure 7: Ergodic Distributions for Selected Years (Phase 3)', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig7, fullfile(figDir, 'fig7_ergodic_selected.png'));
saveas(fig7, fullfile(figDir, 'fig7_ergodic_selected.fig'));

% =====================================================================
%% Print summary statistics
% =====================================================================
fprintf('\n=== Summary Statistics (Phase 3) ===\n\n');

% Find peak recession years
[~, yr2009] = min(abs(yearsList - 2009));
[~, yr1982] = min(abs(yearsList - 1982));
[~, yr2020] = min(abs(yearsList - 2020));
[~, yr1997] = min(abs(yearsList - 1997));
[~, yr2019] = min(abs(yearsList - 2019));

fprintf('Unemployment rates (%%):  1997   2009   2019   2020   1982\n');
for ti = 1:4
    fprintf('  %-7s             %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n', ...
        typeLabels{ti}, Urate(ti,yr1997)*100, Urate(ti,yr2009)*100, ...
        Urate(ti,yr2019)*100, Urate(ti,yr2020)*100, Urate(ti,yr1982)*100);
end
fprintf('  %-7s             %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n', ...
    'Aggreg', agg_Urate(yr1997)*100, agg_Urate(yr2009)*100, ...
    agg_Urate(yr2019)*100, agg_Urate(yr2020)*100, agg_Urate(yr1982)*100);

fprintf('\nEmployment mass:       1997   2009   2019   2020   1982\n');
for ti = 1:4
    fprintf('  %-7s             %5.3f  %5.3f  %5.3f  %5.3f  %5.3f\n', ...
        typeLabels{ti}, E_mass(ti,yr1997), E_mass(ti,yr2009), ...
        E_mass(ti,yr2019), E_mass(ti,yr2020), E_mass(ti,yr1982));
end

fprintf('\nOLF mass:              1997   2009   2019   2020   1982\n');
for ti = 1:4
    fprintf('  %-7s             %5.3f  %5.3f  %5.3f  %5.3f  %5.3f\n', ...
        typeLabels{ti}, N_mass(ti,yr1997), N_mass(ti,yr2009), ...
        N_mass(ti,yr2019), N_mass(ti,yr2020), N_mass(ti,yr1982));
end

fprintf('\nDone. Figures saved to %s\n', figDir);
