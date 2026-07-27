% figures_step3_sequential.m — Compare Sequential vs Existing Phase 3 results
%
% Generates 7 figures from sequential results + overlays with existing Phase 3
% to assess whether sequential warm-starting reduced noise.

clc; clear; close all;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
step1Dir   = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
step3Dir   = fullfile(rewriteDir, 'Step3_FixedOmega');
seqDir     = fullfile(rewriteDir, 'Step3_FixedOmega_Sequential');
figDir     = fullfile(seqDir, 'Figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% Settings
numYears  = 47;
yearsList = setdiff(1976:2024, [1977 1993]);

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
omega_fixed = omf.omega_fixed;

%% Load Phase 1, existing Phase 3, and sequential Phase 3
p1_omega = zeros(numYears, 5);

% Existing Phase 3
ex_ergVecs = zeros(4, 4, numYears);
ex_omega   = zeros(numYears, 5);
ex_fvals   = nan(numYears, 1);

% Sequential Phase 3
sq_ergVecs = zeros(4, 4, numYears);
sq_omega   = zeros(numYears, 5);
sq_fvals   = nan(numYears, 1);
sq_source  = cell(numYears, 1);

for yr = 1:numYears
    % Phase 1
    tmp1 = load(fullfile(step1Dir, sprintf('step1_%d.mat', yr)));
    p1_omega(yr, :) = tmp1.res.omega';

    % Existing Phase 3
    tmp3 = load(fullfile(step3Dir, sprintf('step3_%d.mat', yr)));
    ex_ergVecs(:, :, yr) = tmp3.res.ergVecs;
    ex_omega(yr, :) = tmp3.res.omega';
    ex_fvals(yr) = tmp3.res.SSR2;

    % Sequential Phase 3
    tmpS = load(fullfile(seqDir, sprintf('step3_%d.mat', yr)));
    sq_ergVecs(:, :, yr) = tmpS.res.ergVecs;
    sq_omega(yr, :) = tmpS.res.omega';
    sq_fvals(yr) = tmpS.res.SSR2;
    sq_source{yr} = tmpS.res.source;
end

%% Compute activity masses and unemployment rates for both

% Existing
ex_E = squeeze(ex_ergVecs(3,:,:) + ex_ergVecs(4,:,:));  % 4x47
ex_U = squeeze(ex_ergVecs(2,:,:));
ex_N = squeeze(ex_ergVecs(1,:,:));
ex_Urate = ex_U ./ (ex_U + ex_E);

% Sequential
sq_E = squeeze(sq_ergVecs(3,:,:) + sq_ergVecs(4,:,:));
sq_U = squeeze(sq_ergVecs(2,:,:));
sq_N = squeeze(sq_ergVecs(1,:,:));
sq_Urate = sq_U ./ (sq_U + sq_E);

% Aggregate U-rate
ex_agg_U = zeros(1, numYears);
ex_agg_E = zeros(1, numYears);
sq_agg_U = zeros(1, numYears);
sq_agg_E = zeros(1, numYears);
for yr = 1:numYears
    om_ex = ex_omega(yr,:);
    ex_agg_U(yr) = om_ex(1)*ex_U(1,yr) + om_ex(2)*ex_U(2,yr) + ...
                   om_ex(3)*ex_U(3,yr) + om_ex(4)*ex_U(4,yr);
    ex_agg_E(yr) = om_ex(1)*ex_E(1,yr) + om_ex(2)*ex_E(2,yr) + ...
                   om_ex(3)*ex_E(3,yr) + om_ex(4)*ex_E(4,yr);

    om_sq = sq_omega(yr,:);
    sq_agg_U(yr) = om_sq(1)*sq_U(1,yr) + om_sq(2)*sq_U(2,yr) + ...
                   om_sq(3)*sq_U(3,yr) + om_sq(4)*sq_U(4,yr);
    sq_agg_E(yr) = om_sq(1)*sq_E(1,yr) + om_sq(2)*sq_E(2,yr) + ...
                   om_sq(3)*sq_E(3,yr) + om_sq(4)*sq_E(4,yr);
end
ex_agg_Urate = ex_agg_U ./ (ex_agg_U + ex_agg_E);
sq_agg_Urate = sq_agg_U ./ (sq_agg_U + sq_agg_E);

%% Noise measure: year-to-year absolute changes
fprintf('=== NOISE COMPARISON (mean |year-to-year change|) ===\n');
typeLabels = {'High-N', 'High-U', 'High-E', 'All-E'};
fprintf('%-8s  Existing  Sequential  Reduction\n', 'Type');
fprintf('%-8s  --------  ----------  ---------\n', '');
for ti = 1:4
    ex_noise = mean(abs(diff(ex_Urate(ti,:))));
    sq_noise = mean(abs(diff(sq_Urate(ti,:))));
    fprintf('%-8s  %.4f    %.4f     %.1f%%\n', typeLabels{ti}, ...
        ex_noise, sq_noise, (1-sq_noise/ex_noise)*100);
end
ex_agg_noise = mean(abs(diff(ex_agg_Urate)));
sq_agg_noise = mean(abs(diff(sq_agg_Urate)));
fprintf('%-8s  %.4f    %.4f     %.1f%%\n', 'Aggreg', ...
    ex_agg_noise, sq_agg_noise, (1-sq_agg_noise/ex_agg_noise)*100);

fprintf('\nMean fval: existing=%.4f, sequential=%.4f (delta=%+.4f)\n', ...
    mean(ex_fvals), mean(sq_fvals), mean(sq_fvals)-mean(ex_fvals));

% Helper: add recession bands
    function add_recs(ax, recs)
        yl = ax.YLim;
        for rr = 1:numel(recs)
            rec = recs{rr};
            patch(ax, [rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
                [yl(1), yl(1), yl(2), yl(2)], ...
                [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5, 'HandleVisibility','off');
        end
    end

% =====================================================================
%% FIGURE 1: Type Shares — Phase 1 estimates + Phase 2 trends (same)
% =====================================================================
fig1 = figure('Position', [50 50 900 600], 'Color', 'w');

typeNames = {'High-N', 'High-U', 'High-E', 'All-E', 'All-N'};
trendCols = [3, 4, 5, 6, 2];
typeColors = {cHighN, cHighU, cHighE, cAllE, cAllN};

for i = 1:5
    subplot(2, 3, i);
    hold on;
    yl = [0, max(p1_omega(:,i))*1.3 + 0.01];
    ylim(yl);
    for r = 1:numel(recessions)
        rec = recessions{r};
        patch([rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
              [yl(1), yl(1), yl(2), yl(2)], ...
              [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5, 'HandleVisibility','off');
    end
    plot(yearsList, p1_omega(:,i), 'o', 'Color', typeColors{i}, ...
         'MarkerSize', 4, 'MarkerFaceColor', typeColors{i});
    plot(yearsList, omega_fixed(:, trendCols(i)), '-', 'Color', 'k', ...
         'LineWidth', 2);
    title(typeNames{i}, 'FontSize', 12);
    xlim([1975 2025]);
    if i <= 3, ylabel('Population share'); end
    grid on; box on;
end
sgtitle('Figure 1: Type Shares -- Phase 1 Estimates and Phase 2 Trends', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig1, fullfile(figDir, 'fig1_omega_trends.png'));

% =====================================================================
%% FIGURE 2: Fval comparison (bar chart)
% =====================================================================
fig2 = figure('Position', [50 50 900 400], 'Color', 'w');
hold on;
bar(yearsList, [ex_fvals, sq_fvals - ex_fvals], 'stacked');
colormap([0.7 0.7 0.7; 0.2 0.6 0.2]);
xlabel('Year'); ylabel('fval');
legend('Existing', 'Sequential improvement', 'Location', 'northwest');
title('Figure 2: Fval by Year (Existing vs Sequential)', 'FontSize', 14);
xlim([1975 2025]); grid on; box on;
saveas(fig2, fullfile(figDir, 'fig2_fval_comparison.png'));

% =====================================================================
%% FIGURE 3: Unemployment Rate — THE KEY COMPARISON
%  Existing = dashed thin, Sequential = solid thick
% =====================================================================
fig3 = figure('Position', [50 50 1000 550], 'Color', 'w');
hold on;

yl = [0 min(100, max([ex_Urate(:); sq_Urate(:)])*110)];
ylim(yl);
for r = 1:numel(recessions)
    rec = recessions{r};
    patch([rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
          [yl(1), yl(1), yl(2), yl(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5, 'HandleVisibility','off');
end

% Existing (dashed, thin)
plot(yearsList, ex_Urate(1,:)*100, '--', 'Color', cHighN, 'LineWidth', 1);
plot(yearsList, ex_Urate(2,:)*100, '--', 'Color', cHighU, 'LineWidth', 1);
plot(yearsList, ex_Urate(3,:)*100, '--', 'Color', cHighE, 'LineWidth', 1);
plot(yearsList, ex_Urate(4,:)*100, '--', 'Color', cAllE, 'LineWidth', 1);
plot(yearsList, ex_agg_Urate*100,  '--k', 'LineWidth', 1);

% Sequential (solid, thick)
h1 = plot(yearsList, sq_Urate(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 2, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
h2 = plot(yearsList, sq_Urate(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 2, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
h3 = plot(yearsList, sq_Urate(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 2, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
h4 = plot(yearsList, sq_Urate(4,:)*100, '-v', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');
h5 = plot(yearsList, sq_agg_Urate*100, '-k', 'LineWidth', 2.5, 'DisplayName', 'Model aggregate');

% Dummy for legend
h6 = plot(NaN, NaN, '--k', 'LineWidth', 1, 'DisplayName', 'Existing (dashed)');

legend([h1 h2 h3 h4 h5 h6], 'Location', 'northeast', 'FontSize', 10);
xlabel('Year', 'FontSize', 12);
ylabel('Unemployment Rate, %', 'FontSize', 12);
title('Figure 3: U-Rate by Type -- Sequential (solid) vs Existing (dashed)', 'FontSize', 13);
xlim([1975 2025]); grid on; box on;
set(gca, 'FontSize', 11);
saveas(fig3, fullfile(figDir, 'fig3_urate_comparison.png'));

% =====================================================================
%% FIGURE 4: U-Rate Mover Types Only — Side by side
% =====================================================================
fig4 = figure('Position', [50 50 1100 450], 'Color', 'w');

% Panel A: Existing
subplot(1,2,1); hold on;
yl4 = [0 max([ex_Urate(1:3,:)*100, sq_Urate(1:3,:)*100], [], 'all')*1.1];
ylim(yl4);
for r = 1:numel(recessions)
    rec = recessions{r};
    patch([rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
          [yl4(1), yl4(1), yl4(2), yl4(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5, 'HandleVisibility','off');
end
plot(yearsList, ex_Urate(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, ex_Urate(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, ex_Urate(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
legend('Location', 'northeast', 'FontSize', 10);
xlabel('Year'); ylabel('Unemployment Rate, %');
title('(A) Existing Phase 3', 'FontSize', 12);
xlim([1975 2025]); grid on; box on;

% Panel B: Sequential
subplot(1,2,2); hold on;
ylim(yl4);
for r = 1:numel(recessions)
    rec = recessions{r};
    patch([rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
          [yl4(1), yl4(1), yl4(2), yl4(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5, 'HandleVisibility','off');
end
plot(yearsList, sq_Urate(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, sq_Urate(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, sq_Urate(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
legend('Location', 'northeast', 'FontSize', 10);
xlabel('Year'); ylabel('Unemployment Rate, %');
title('(B) Sequential Phase 3', 'FontSize', 12);
xlim([1975 2025]); grid on; box on;

sgtitle('Figure 4: U-Rate Mover Types -- Existing vs Sequential', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig4, fullfile(figDir, 'fig4_urate_movers_sidebyside.png'));

% =====================================================================
%% FIGURE 5: Ergodic E-mass by type (sequential only)
% =====================================================================
fig5 = figure('Position', [50 50 800 450], 'Color', 'w');
hold on;
yl = [0 1.05]; ylim(yl);
for r = 1:numel(recessions)
    rec = recessions{r};
    patch([rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
          [yl(1), yl(1), yl(2), yl(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5, 'HandleVisibility','off');
end
plot(yearsList, sq_E(1,:), '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, sq_E(2,:), '-o', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, sq_E(3,:), '-o', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
plot(yearsList, sq_E(4,:), '-o', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');
legend('Location', 'southwest', 'FontSize', 10);
xlabel('Year'); ylabel('Ergodic E-mass (STJ + LTJ)');
title('Figure 5: Employment Mass by Type (Sequential)', 'FontSize', 14);
xlim([1975 2025]); grid on; box on;
saveas(fig5, fullfile(figDir, 'fig5_emass_bytype.png'));

% =====================================================================
%% FIGURE 6: Ergodic U-mass by type (sequential only)
% =====================================================================
fig6 = figure('Position', [50 50 800 450], 'Color', 'w');
hold on;
yl = [0 max(sq_U(:))*1.15]; ylim(yl);
for r = 1:numel(recessions)
    rec = recessions{r};
    patch([rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
          [yl(1), yl(1), yl(2), yl(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5, 'HandleVisibility','off');
end
plot(yearsList, sq_U(1,:), '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, sq_U(2,:), '-o', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, sq_U(3,:), '-o', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
plot(yearsList, sq_U(4,:), '-o', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');
legend('Location', 'northeast', 'FontSize', 10);
xlabel('Year'); ylabel('Ergodic U-mass');
title('Figure 6: Unemployment Mass by Type (Sequential)', 'FontSize', 14);
xlim([1975 2025]); grid on; box on;
saveas(fig6, fullfile(figDir, 'fig6_umass_bytype.png'));

% =====================================================================
%% FIGURE 7: Ergodic Distributions — Stacked bars for selected years
% =====================================================================
selectedYears = [1978, 1982, 1997, 2007, 2009, 2019, 2020, 2024];
selIdx = arrayfun(@(y) find(yearsList == y), selectedYears);

fig7 = figure('Position', [50 50 1000 400], 'Color', 'w');
actLabels = {'OLF', 'U', 'E'};

for ti = 1:4
    subplot(1, 4, ti);
    barData = zeros(numel(selectedYears), 3);
    for si = 1:numel(selectedYears)
        barData(si, :) = [sq_N(ti, selIdx(si)), ...
                          sq_U(ti, selIdx(si)), ...
                          sq_E(ti, selIdx(si))];
    end
    b = bar(barData, 'stacked');
    b(1).FaceColor = cHighN;
    b(2).FaceColor = cHighU;
    b(3).FaceColor = cHighE;
    set(gca, 'XTickLabel', selectedYears, 'XTickLabelRotation', 45);
    ylim([0 1.05]);
    title(typeLabels{ti}, 'FontSize', 12);
    if ti == 1
        ylabel('Ergodic share');
        legend(actLabels, 'Location', 'northwest', 'FontSize', 8);
    end
    grid on; box on;
end
sgtitle('Figure 7: Ergodic Distributions for Selected Years (Sequential)', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig7, fullfile(figDir, 'fig7_ergodic_selected.png'));

% =====================================================================
%% Print summary statistics
% =====================================================================
fprintf('\n=== Summary Statistics (Sequential Phase 3) ===\n\n');

[~, yr2009] = min(abs(yearsList - 2009));
[~, yr1982] = min(abs(yearsList - 1982));
[~, yr2020] = min(abs(yearsList - 2020));
[~, yr1997] = min(abs(yearsList - 1997));
[~, yr2019] = min(abs(yearsList - 2019));

fprintf('Unemployment rates (%%):  1997   2009   2019   2020   1982\n');
for ti = 1:4
    fprintf('  %-7s             %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n', ...
        typeLabels{ti}, sq_Urate(ti,yr1997)*100, sq_Urate(ti,yr2009)*100, ...
        sq_Urate(ti,yr2019)*100, sq_Urate(ti,yr2020)*100, sq_Urate(ti,yr1982)*100);
end
fprintf('  %-7s             %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n', ...
    'Aggreg', sq_agg_Urate(yr1997)*100, sq_agg_Urate(yr2009)*100, ...
    sq_agg_Urate(yr2019)*100, sq_agg_Urate(yr2020)*100, sq_agg_Urate(yr1982)*100);

fprintf('\nSource breakdown:\n');
sources_unique = unique(sq_source);
for s = 1:numel(sources_unique)
    n = sum(strcmp(sq_source, sources_unique{s}));
    fprintf('  %s: %d years\n', sources_unique{s}, n);
end

fprintf('\nImproved years (delta < -1e-4): %d/%d\n', ...
    sum(sq_fvals < ex_fvals - 1e-4), numYears);

fprintf('\nDone. Figures saved to %s\n', figDir);
