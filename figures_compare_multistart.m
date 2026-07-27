% figures_compare_multistart.m — Compare original vs multi-start Step 1 results
%
% Produces side-by-side comparison figures:
%   Fig C1: Type shares (omegas) — old vs new, 5 panels
%   Fig C2: Unemployment rate by mover type — old vs new overlay
%   Fig C3: Ergodic distributions — old vs new, 4 panels x 3 types
%   Fig C4: Fval comparison across years
%   Fig C5: New results only — all type shares (publication figure)
%   Fig C6: New results only — unemployment rate by mover type
%   Fig C7: New results only — ergodic distributions (4-panel)
%   Fig C8: New results only — unemployment rate by type + aggregate

rewriteDir = fileparts(mfilename('fullpath'));
oldDir = fullfile(rewriteDir, 'Step1_FreeOmega');
newDir = fullfile(rewriteDir, 'Step1_FreeOmega_Improving');
figDir = fullfile(rewriteDir, 'Step1_FreeOmega_Improving', 'Figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

yearsList = setdiff(1976:2024, [1977 1993]);
numYears = numel(yearsList);
Theta = 3; ThetaTot = 5;

%% NBER recession periods
recessions = [1980 1980; 1981 1982; 1990 1991; 2001 2001; 2008 2009; 2020 2020];

%% Load both sets of results
omega_old = zeros(ThetaTot, numYears);
omega_new = zeros(ThetaTot, numYears);
erg_old   = zeros(4, Theta, numYears);
erg_new   = zeros(4, Theta, numYears);
ergAgg_old = zeros(4, numYears);
ergAgg_new = zeros(4, numYears);
fval_old  = zeros(1, numYears);
fval_new  = zeros(1, numYears);

for yr = 1:numYears
    % Old results
    fname = fullfile(oldDir, sprintf('step1_%d.mat', yr));
    if isfile(fname)
        tmp = load(fname); r = tmp.res;
        omega_old(:,yr) = r.omega;
        erg_old(:,:,yr) = r.ergVecs;
        ergAgg_old(:,yr) = r.ergVecsAll;
        fval_old(yr) = r.SSR2;
    end
    % New results
    fname = fullfile(newDir, sprintf('step1_%d.mat', yr));
    if isfile(fname)
        tmp = load(fname); r = tmp.res;
        omega_new(:,yr) = r.omega;
        erg_new(:,:,yr) = r.ergVecs;
        ergAgg_new(:,yr) = r.ergVecsAll;
        fval_new(yr) = r.SSR2;
    end
end

%% Color scheme
c_HighN = [0.00 0.45 0.74];   % blue
c_HighU = [0.85 0.33 0.10];   % red-orange
c_HighE = [0.47 0.67 0.19];   % green
c_AllE  = [0.49 0.18 0.56];   % purple
c_AllN  = [0.64 0.08 0.18];   % dark red
c_Agg   = [0.30 0.30 0.30];   % dark gray

typeIdx   = [4, 3, 2, 1, 5];  % omega order: All-E, High-E, High-U, High-N, All-N
typeNames = {'All-E','High-E','High-U','High-N','All-N'};
typeCol   = {c_AllE, c_HighE, c_HighU, c_HighN, c_AllN};

%% ========================================================================
%  FIG C1: Type shares — old (dashed) vs new (solid), 5 panels
%  ========================================================================
figC1 = figure('Color','w','Position',[50 50 1400 900]);
for k = 1:5
    ax = subplot(3,2,k); hold(ax,'on');
    idx = typeIdx(k);
    ym = max(max(omega_old(idx,:)), max(omega_new(idx,:)))*1.15 + 0.02;
    add_recession_bars(ax, recessions, [0 ym]);
    plot(ax, yearsList, omega_old(idx,:), '--', 'LineWidth', 1.5, 'Color', [0.6 0.6 0.6]);
    plot(ax, yearsList, omega_new(idx,:), '-',  'LineWidth', 2.0, 'Color', typeCol{k});
    title(ax, typeNames{k}, 'FontSize', 12);
    xlim(ax, [1976 2024]); ylim(ax, [0 ym]);
    set(ax, 'FontSize', 10, 'Box', 'on', 'YGrid', 'on');
    if k == 1
        legend(ax, {'Original','Multi-start'}, 'Location','best', 'FontSize', 9);
    end
    if k >= 4, xlabel(ax, 'Year'); end
    ylabel(ax, '\omega');
end
% Use the 6th subplot for annotation
ax6 = subplot(3,2,6); axis(ax6, 'off');
improved = find(fval_new < fval_old - 1e-6);
text(ax6, 0.05, 0.9, sprintf('Years with improved fval: %d / %d', numel(improved), numYears), ...
    'FontSize', 12, 'FontWeight', 'bold');
impText = '';
for j = 1:numel(improved)
    impText = [impText, sprintf('%d (%.4f), ', yearsList(improved(j)), fval_old(improved(j))-fval_new(improved(j)))];
    if mod(j,4) == 0, impText = [impText, newline]; end
end
text(ax6, 0.05, 0.6, impText, 'FontSize', 9, 'VerticalAlignment', 'top');
sgtitle(figC1, 'Type Shares: Original (dashed gray) vs Multi-start (solid color)', 'FontSize', 14);
exportgraphics(figC1, fullfile(figDir, 'FigC1_omegas_comparison.png'), 'Resolution', 300);

%% ========================================================================
%  FIG C2: Unemployment rate by mover type — old vs new
%  ========================================================================
urate_old = zeros(Theta, numYears);
urate_new = zeros(Theta, numYears);
for yr = 1:numYears
    for th = 1:Theta
        ev = erg_old(:,th,yr);
        d = ev(2)+ev(3)+ev(4);
        if d > 0, urate_old(th,yr) = ev(2)/d; end
        ev = erg_new(:,th,yr);
        d = ev(2)+ev(3)+ev(4);
        if d > 0, urate_new(th,yr) = ev(2)/d; end
    end
end

moverNames = {'High-N','High-U','High-E'};
moverCol   = {c_HighN, c_HighU, c_HighE};

figC2 = figure('Color','w','Position',[100 100 1200 500]);
for th = 1:Theta
    ax = subplot(1,3,th); hold(ax,'on');
    add_recession_bars(ax, recessions, [0 1]);
    plot(ax, yearsList, urate_old(th,:), '--', 'LineWidth', 1.5, 'Color', [0.6 0.6 0.6]);
    plot(ax, yearsList, urate_new(th,:), '-',  'LineWidth', 2.0, 'Color', moverCol{th});
    title(ax, moverNames{th}, 'FontSize', 12);
    xlim(ax, [1976 2024]); ylim(ax, [0 1]);
    set(ax, 'FontSize', 10, 'Box', 'on', 'YGrid', 'on');
    xlabel(ax, 'Year'); ylabel(ax, 'U-rate');
    if th == 1
        legend(ax, {'Original','Multi-start'}, 'Location','best', 'FontSize', 9);
    end
end
sgtitle(figC2, 'Unemployment Rate by Mover Type: Original vs Multi-start', 'FontSize', 14);
exportgraphics(figC2, fullfile(figDir, 'FigC2_urate_comparison.png'), 'Resolution', 300);

%% ========================================================================
%  FIG C3: Ergodic distributions — old vs new, 4 states x 3 types
%  ========================================================================
stateLabels = {'N (OLF)', 'U', 'E short-term', 'E long-term'};

figC3 = figure('Color','w','Position',[50 50 1400 900]);
for s = 1:4
    for th = 1:Theta
        pIdx = (s-1)*3 + th;
        ax = subplot(4,3,pIdx); hold(ax,'on');
        add_recession_bars(ax, recessions, [0 1]);
        plot(ax, yearsList, squeeze(erg_old(s,th,:))', '--', 'LineWidth', 1.3, 'Color', [0.6 0.6 0.6]);
        plot(ax, yearsList, squeeze(erg_new(s,th,:))', '-',  'LineWidth', 1.8, 'Color', moverCol{th});
        if s == 1, title(ax, moverNames{th}, 'FontSize', 11); end
        if th == 1, ylabel(ax, stateLabels{s}, 'FontSize', 10); end
        xlim(ax, [1976 2024]);
        set(ax, 'FontSize', 8, 'Box', 'on', 'YGrid', 'on');
        if s == 4, xlabel(ax, 'Year'); end
    end
end
sgtitle(figC3, 'Ergodic Distributions: Original (dashed) vs Multi-start (solid)', 'FontSize', 14);
exportgraphics(figC3, fullfile(figDir, 'FigC3_ergodic_comparison.png'), 'Resolution', 300);

%% ========================================================================
%  FIG C4: Fval comparison
%  ========================================================================
figC4 = figure('Color','w','Position',[100 100 900 400]);
ax = axes(figC4); hold(ax,'on');
delta = fval_old - fval_new;  % positive = new is better
bar(ax, yearsList, delta, 0.7, 'FaceColor', [0.2 0.6 0.2], 'EdgeColor', 'none');
% Mark years that got worse
worse = delta < -1e-6;
if any(worse)
    bar(ax, yearsList(worse), delta(worse), 0.7, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none');
end
yline(ax, 0, 'k-', 'LineWidth', 0.5);
xlabel(ax, 'Year', 'FontSize', 12);
ylabel(ax, '\Delta fval (old - new)', 'FontSize', 12);
title(ax, 'Improvement in objective: positive = multi-start found better solution', 'FontSize', 13);
set(ax, 'FontSize', 10, 'Box', 'on', 'YGrid', 'on');
xlim(ax, [1975 2025]);
exportgraphics(figC4, fullfile(figDir, 'FigC4_fval_comparison.png'), 'Resolution', 300);

%% ========================================================================
%  FIG C5: New results — all type shares (publication figure)
%  ========================================================================
figC5 = figure('Color','w','Position',[100 100 900 450]);
ax = axes(figC5); hold(ax,'on');
add_recession_bars(ax, recessions, [0 1]);
h = gobjects(5,1);
for k = 1:5
    h(k) = plot(ax, yearsList, omega_new(typeIdx(k),:), '-', 'LineWidth', 2.2, 'Color', typeCol{k});
end
legend(h, typeNames, 'Location','east', 'FontSize',11);
xlabel(ax, 'Year', 'FontSize',12);
ylabel(ax, 'Type share (\omega)', 'FontSize',12);
title(ax, 'Step 1: Estimated type shares (multi-start)', 'FontSize',14);
xlim(ax, [1976 2024]); ylim(ax, [0 1]);
set(ax, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(figC5, fullfile(figDir, 'FigC5_omegas_new.png'), 'Resolution', 300);

%% ========================================================================
%  FIG C6: New results — unemployment rate by mover type
%  ========================================================================
figC6 = figure('Color','w','Position',[100 100 900 450]);
ax = axes(figC6); hold(ax,'on');
add_recession_bars(ax, recessions, [0 1]);
h6 = gobjects(3,1);
for th = 1:Theta
    h6(th) = plot(ax, yearsList, urate_new(th,:), '-', 'LineWidth', 2.2, 'Color', moverCol{th});
end
legend(h6, moverNames, 'Location','northeast', 'FontSize',11);
xlabel(ax, 'Year', 'FontSize',12);
ylabel(ax, 'Unemployment rate', 'FontSize',12);
title(ax, 'Step 1: Unemployment rate by mover type (multi-start)', 'FontSize',14);
xlim(ax, [1976 2024]);
set(ax, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(figC6, fullfile(figDir, 'FigC6_urate_mover_new.png'), 'Resolution', 300);

%% ========================================================================
%  FIG C7: New results — ergodic distributions (4-panel)
%  ========================================================================
figC7 = figure('Color','w','Position',[100 100 1200 800]);
for s = 1:4
    ax = subplot(2,2,s); hold(ax,'on');
    add_recession_bars(ax, recessions, [0 1]);
    h7 = gobjects(Theta,1);
    for th = 1:Theta
        h7(th) = plot(ax, yearsList, squeeze(erg_new(s,th,:))', '-', ...
            'LineWidth', 2, 'Color', moverCol{th});
    end
    if s == 1, legend(h7, moverNames, 'Location','best', 'FontSize',9); end
    title(ax, sprintf('Ergodic share: %s', stateLabels{s}), 'FontSize',12);
    xlabel(ax, 'Year'); ylabel(ax, 'Share');
    xlim(ax, [1976 2024]);
    set(ax, 'FontSize',10, 'Box','on', 'YGrid','on');
end
sgtitle(figC7, 'Step 1: Ergodic distributions by type (multi-start)', 'FontSize',14);
exportgraphics(figC7, fullfile(figDir, 'FigC7_ergodic_new.png'), 'Resolution', 300);

%% ========================================================================
%  FIG C8: New results — unemployment rate by type + aggregate
%  ========================================================================
urate_agg_new = zeros(1, numYears);
for yr = 1:numYears
    ea = ergAgg_new(:,yr);
    d = ea(2)+ea(3)+ea(4);
    if d > 0, urate_agg_new(yr) = ea(2)/d; end
end

figC8 = figure('Color','w','Position',[100 100 900 450]);
ax = axes(figC8); hold(ax,'on');
add_recession_bars(ax, recessions, [0 1]);
h8 = gobjects(4,1);
for th = 1:Theta
    h8(th) = plot(ax, yearsList, urate_new(th,:), '-', 'LineWidth', 1.8, 'Color', moverCol{th});
end
h8(4) = plot(ax, yearsList, urate_agg_new, '-', 'LineWidth', 2.5, 'Color', c_Agg);
legend(h8, [moverNames, {'Aggregate'}], 'Location','northeast', 'FontSize',11);
xlabel(ax, 'Year', 'FontSize',12);
ylabel(ax, 'Unemployment rate', 'FontSize',12);
title(ax, 'Step 1: Unemployment rate by type + aggregate (multi-start)', 'FontSize',14);
xlim(ax, [1976 2024]);
set(ax, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(figC8, fullfile(figDir, 'FigC8_urate_all_new.png'), 'Resolution', 300);

%% Print summary
fprintf('\n=== COMPARISON SUMMARY ===\n');
fprintf('Years with improved fval: %d / %d\n', sum(fval_new < fval_old - 1e-6), numYears);
fprintf('Years with worse fval:    %d / %d\n', sum(fval_new > fval_old + 1e-6), numYears);
fprintf('Max improvement: %.4f (year %d)\n', max(delta), yearsList(delta == max(delta)));
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
