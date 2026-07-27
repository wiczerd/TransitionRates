% build_best_of_both.m — Combine original and multi-start results
%
% For each year, picks the solution with the lower fval (SSR2).
% Saves combined results to Step1_FreeOmega_Best/
% Then generates comparison figures across all three sets.

clc; clear;

rewriteDir = fileparts(mfilename('fullpath'));
oldDir  = fullfile(rewriteDir, 'Step1_FreeOmega');
msDir   = fullfile(rewriteDir, 'Step1_FreeOmega_Improving');
bestDir = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
figDir  = fullfile(bestDir, 'Figures');
if ~exist(bestDir, 'dir'), mkdir(bestDir); end
if ~exist(figDir, 'dir'), mkdir(figDir); end

yearsList = setdiff(1976:2024, [1977 1993]);
numYears  = numel(yearsList);
Theta = 3; ThetaTot = 5;

%% NBER recession periods
recessions = [1980 1980; 1981 1982; 1990 1991; 2001 2001; 2008 2009; 2020 2020];

%% Compare and pick best for each year
fval_old  = nan(1, numYears);
fval_ms   = nan(1, numYears);
fval_best = nan(1, numYears);
source    = cell(1, numYears);  % 'old' or 'ms'

omega_old  = zeros(ThetaTot, numYears);
omega_ms   = zeros(ThetaTot, numYears);
omega_best = zeros(ThetaTot, numYears);

erg_old    = zeros(4, Theta, numYears);
erg_ms     = zeros(4, Theta, numYears);
erg_best   = zeros(4, Theta, numYears);

ergAgg_old  = zeros(4, numYears);
ergAgg_ms   = zeros(4, numYears);
ergAgg_best = zeros(4, numYears);

fprintf('=== Building best-of-both results ===\n\n');
fprintf('%-6s %-10s %-10s %-10s %-8s\n', 'Year', 'Old fval', 'MS fval', 'Best fval', 'Source');
fprintf('%s\n', repmat('-', 1, 50));

for yr = 1:numYears
    % Load old
    tmp = load(fullfile(oldDir, sprintf('step1_%d.mat', yr)));
    r_old = tmp.res;
    fval_old(yr) = r_old.SSR2;
    omega_old(:,yr) = r_old.omega;
    erg_old(:,:,yr) = r_old.ergVecs;
    ergAgg_old(:,yr) = r_old.ergVecsAll;

    % Load multi-start
    tmp = load(fullfile(msDir, sprintf('step1_%d.mat', yr)));
    r_ms = tmp.res;
    fval_ms(yr) = r_ms.SSR2;
    omega_ms(:,yr) = r_ms.omega;
    erg_ms(:,:,yr) = r_ms.ergVecs;
    ergAgg_ms(:,yr) = r_ms.ergVecsAll;

    % Pick best
    if fval_old(yr) <= fval_ms(yr)
        source{yr} = 'old';
        best_res = r_old;
    else
        source{yr} = 'ms';
        best_res = r_ms;
    end

    fval_best(yr) = best_res.SSR2;
    omega_best(:,yr) = best_res.omega;
    erg_best(:,:,yr) = best_res.ergVecs;
    ergAgg_best(:,yr) = best_res.ergVecsAll;

    % Save best to new folder
    res = best_res;
    save(fullfile(bestDir, sprintf('step1_%d.mat', yr)), 'res');

    flag = '';
    if ~strcmp(source{yr}, 'old')
        flag = ' <-- MS wins';
    end
    fprintf('%-6d %-10.4f %-10.4f %-10.4f %-8s%s\n', ...
        yearsList(yr), fval_old(yr), fval_ms(yr), fval_best(yr), source{yr}, flag);
end

nFromOld = sum(strcmp(source, 'old'));
nFromMS  = sum(strcmp(source, 'ms'));
fprintf('\n=== Summary ===\n');
fprintf('From original:    %d / %d years\n', nFromOld, numYears);
fprintf('From multi-start: %d / %d years\n', nFromMS, numYears);
fprintf('Total fval improvement: %.4f\n', sum(fval_old) - sum(fval_best));

%% Save summary table
T = table(yearsList(:), fval_old(:), fval_ms(:), fval_best(:), source(:), ...
    'VariableNames', {'Year','Fval_Old','Fval_MS','Fval_Best','Source'});
writetable(T, fullfile(bestDir, 'best_of_both_summary.csv'));

%% ========================================================================
%  Color scheme
%  ========================================================================
c_HighN = [0.00 0.45 0.74];   % blue
c_HighU = [0.85 0.33 0.10];   % red-orange
c_HighE = [0.47 0.67 0.19];   % green
c_AllE  = [0.49 0.18 0.56];   % purple
c_AllN  = [0.64 0.08 0.18];   % dark red
c_Agg   = [0.30 0.30 0.30];   % dark gray

typeIdx   = [4, 3, 2, 1, 5];  % omega order: All-E, High-E, High-U, High-N, All-N
typeNames = {'All-E','High-E','High-U','High-N','All-N'};
typeCol   = {c_AllE, c_HighE, c_HighU, c_HighN, c_AllN};
moverNames = {'High-N','High-U','High-E'};
moverCol   = {c_HighN, c_HighU, c_HighE};

%% Compute u-rates for all three sets
urate_old  = zeros(Theta, numYears);
urate_ms   = zeros(Theta, numYears);
urate_best = zeros(Theta, numYears);
for yr = 1:numYears
    for th = 1:Theta
        ev = erg_old(:,th,yr);  d = ev(2)+ev(3)+ev(4);
        if d > 0, urate_old(th,yr) = ev(2)/d; end
        ev = erg_ms(:,th,yr);   d = ev(2)+ev(3)+ev(4);
        if d > 0, urate_ms(th,yr) = ev(2)/d; end
        ev = erg_best(:,th,yr); d = ev(2)+ev(3)+ev(4);
        if d > 0, urate_best(th,yr) = ev(2)/d; end
    end
end

%% ========================================================================
%  FIG B1: Type shares — old (dotted gray), MS (dashed gray), best (solid)
%  ========================================================================
figB1 = figure('Color','w','Position',[50 50 1400 900]);
for k = 1:5
    ax = subplot(3,2,k); hold(ax,'on');
    idx = typeIdx(k);
    ym = max([max(omega_old(idx,:)), max(omega_ms(idx,:)), max(omega_best(idx,:))])*1.15 + 0.02;
    add_recession_bars(ax, recessions, [0 ym]);
    plot(ax, yearsList, omega_old(idx,:),  ':', 'LineWidth', 1.2, 'Color', [0.7 0.7 0.7]);
    plot(ax, yearsList, omega_ms(idx,:),  '--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
    plot(ax, yearsList, omega_best(idx,:), '-', 'LineWidth', 2.2, 'Color', typeCol{k});
    title(ax, typeNames{k}, 'FontSize', 12);
    xlim(ax, [1976 2024]); ylim(ax, [0 ym]);
    set(ax, 'FontSize', 10, 'Box', 'on', 'YGrid', 'on');
    if k == 1
        legend(ax, {'Original','Multi-start','Best'}, 'Location','best', 'FontSize', 9);
    end
    if k >= 4, xlabel(ax, 'Year'); end
    ylabel(ax, '\omega');
end
ax6 = subplot(3,2,6); axis(ax6, 'off');
text(ax6, 0.05, 0.9, sprintf('Best-of-both: %d old + %d multi-start', nFromOld, nFromMS), ...
    'FontSize', 13, 'FontWeight', 'bold');
msYears = yearsList(strcmp(source, 'ms'));
text(ax6, 0.05, 0.7, sprintf('MS wins: %s', mat2str(msYears)), 'FontSize', 10);
sgtitle(figB1, 'Type Shares: Original (dotted) / Multi-start (dashed) / Best (solid)', 'FontSize', 14);
exportgraphics(figB1, fullfile(figDir, 'FigB1_omegas_3way.png'), 'Resolution', 300);

%% ========================================================================
%  FIG B2: U-rate by mover type — 3 panels, 3 sets each
%  ========================================================================
figB2 = figure('Color','w','Position',[100 100 1200 500]);
for th = 1:Theta
    ax = subplot(1,3,th); hold(ax,'on');
    add_recession_bars(ax, recessions, [0 1]);
    plot(ax, yearsList, urate_old(th,:),  ':', 'LineWidth', 1.2, 'Color', [0.7 0.7 0.7]);
    plot(ax, yearsList, urate_ms(th,:),  '--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
    plot(ax, yearsList, urate_best(th,:), '-', 'LineWidth', 2.0, 'Color', moverCol{th});
    title(ax, moverNames{th}, 'FontSize', 12);
    xlim(ax, [1976 2024]); ylim(ax, [0 1]);
    set(ax, 'FontSize', 10, 'Box', 'on', 'YGrid', 'on');
    xlabel(ax, 'Year'); ylabel(ax, 'U-rate');
    if th == 1
        legend(ax, {'Original','Multi-start','Best'}, 'Location','best', 'FontSize', 9);
    end
end
sgtitle(figB2, 'Unemployment Rate: Original / Multi-start / Best', 'FontSize', 14);
exportgraphics(figB2, fullfile(figDir, 'FigB2_urate_3way.png'), 'Resolution', 300);

%% ========================================================================
%  FIG B3: Ergodic distributions — 4 states x 3 types, best only
%  ========================================================================
stateLabels = {'N (OLF)', 'U', 'E short-term', 'E long-term'};
figB3 = figure('Color','w','Position',[50 50 1200 800]);
for s = 1:4
    ax = subplot(2,2,s); hold(ax,'on');
    add_recession_bars(ax, recessions, [0 1]);
    h = gobjects(Theta,1);
    for th = 1:Theta
        h(th) = plot(ax, yearsList, squeeze(erg_best(s,th,:))', '-', ...
            'LineWidth', 2, 'Color', moverCol{th});
    end
    if s == 1, legend(h, moverNames, 'Location','best', 'FontSize',9); end
    title(ax, sprintf('Ergodic share: %s', stateLabels{s}), 'FontSize',12);
    xlabel(ax, 'Year'); ylabel(ax, 'Share');
    xlim(ax, [1976 2024]);
    set(ax, 'FontSize',10, 'Box','on', 'YGrid','on');
end
sgtitle(figB3, 'Step 1: Ergodic distributions by type (best-of-both)', 'FontSize',14);
exportgraphics(figB3, fullfile(figDir, 'FigB3_ergodic_best.png'), 'Resolution', 300);

%% ========================================================================
%  FIG B4: Best-of-both — all type shares (publication figure)
%  ========================================================================
figB4 = figure('Color','w','Position',[100 100 900 450]);
ax = axes(figB4); hold(ax,'on');
add_recession_bars(ax, recessions, [0 1]);
h = gobjects(5,1);
for k = 1:5
    h(k) = plot(ax, yearsList, omega_best(typeIdx(k),:), '-', 'LineWidth', 2.2, 'Color', typeCol{k});
end
legend(h, typeNames, 'Location','east', 'FontSize',11);
xlabel(ax, 'Year', 'FontSize',12);
ylabel(ax, 'Type share (\omega)', 'FontSize',12);
title(ax, 'Step 1: Estimated type shares (best-of-both)', 'FontSize',14);
xlim(ax, [1976 2024]); ylim(ax, [0 1]);
set(ax, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(figB4, fullfile(figDir, 'FigB4_omegas_best.png'), 'Resolution', 300);

%% ========================================================================
%  FIG B5: Best-of-both — u-rate by mover type
%  ========================================================================
figB5 = figure('Color','w','Position',[100 100 900 450]);
ax = axes(figB5); hold(ax,'on');
add_recession_bars(ax, recessions, [0 1]);
h5 = gobjects(3,1);
for th = 1:Theta
    h5(th) = plot(ax, yearsList, urate_best(th,:), '-', 'LineWidth', 2.2, 'Color', moverCol{th});
end
legend(h5, moverNames, 'Location','northeast', 'FontSize',11);
xlabel(ax, 'Year', 'FontSize',12);
ylabel(ax, 'Unemployment rate', 'FontSize',12);
title(ax, 'Step 1: Unemployment rate by mover type (best-of-both)', 'FontSize',14);
xlim(ax, [1976 2024]);
set(ax, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(figB5, fullfile(figDir, 'FigB5_urate_mover_best.png'), 'Resolution', 300);

%% ========================================================================
%  FIG B6: Best-of-both — u-rate by type + aggregate
%  ========================================================================
urate_agg = zeros(1, numYears);
for yr = 1:numYears
    ea = ergAgg_best(:,yr);
    d = ea(2)+ea(3)+ea(4);
    if d > 0, urate_agg(yr) = ea(2)/d; end
end

figB6 = figure('Color','w','Position',[100 100 900 450]);
ax = axes(figB6); hold(ax,'on');
add_recession_bars(ax, recessions, [0 1]);
h6 = gobjects(4,1);
for th = 1:Theta
    h6(th) = plot(ax, yearsList, urate_best(th,:), '-', 'LineWidth', 1.8, 'Color', moverCol{th});
end
h6(4) = plot(ax, yearsList, urate_agg, '-', 'LineWidth', 2.5, 'Color', c_Agg);
legend(h6, [moverNames, {'Aggregate'}], 'Location','northeast', 'FontSize',11);
xlabel(ax, 'Year', 'FontSize',12);
ylabel(ax, 'Unemployment rate', 'FontSize',12);
title(ax, 'Step 1: Unemployment rate by type + aggregate (best-of-both)', 'FontSize',14);
xlim(ax, [1976 2024]);
set(ax, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(figB6, fullfile(figDir, 'FigB6_urate_all_best.png'), 'Resolution', 300);

%% ========================================================================
%  FIG B7: Fval comparison bar chart
%  ========================================================================
figB7 = figure('Color','w','Position',[100 100 900 400]);
ax = axes(figB7); hold(ax,'on');
delta = fval_old - fval_best;  % always >= 0 by construction
bar(ax, yearsList, delta, 0.7, 'FaceColor', [0.2 0.6 0.2], 'EdgeColor', 'none');
xlabel(ax, 'Year', 'FontSize', 12);
ylabel(ax, '\Delta fval (old - best)', 'FontSize', 12);
title(ax, sprintf('Best-of-both improvement over original (%d years improved)', nFromMS), 'FontSize', 13);
set(ax, 'FontSize', 10, 'Box', 'on', 'YGrid', 'on');
xlim(ax, [1975 2025]);
exportgraphics(figB7, fullfile(figDir, 'FigB7_fval_improvement.png'), 'Resolution', 300);

fprintf('\nAll figures saved to:\n  %s\n', figDir);

%% ========================================================================
%  Helper: recession shading
%  ========================================================================
function add_recession_bars(ax, recessions, ylims)
    for r = 1:size(recessions,1)
        x1 = recessions(r,1) - 0.5;
        x2 = recessions(r,2) + 0.5;
        fill(ax, [x1 x2 x2 x1], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
            [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.5);
    end
    children = get(ax, 'Children');
    nRec = size(recessions,1);
    if numel(children) >= nRec
        set(ax, 'Children', [children(nRec+1:end); children(1:nRec)]);
    end
end
