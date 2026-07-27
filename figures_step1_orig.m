% figures_step1_orig.m — Publication-quality figures for Step 1 (original code output)
%
% Produces:
%   Fig 1: All type shares over time
%   Fig 2a-2e: Individual type shares (5 separate figures)
%   Fig 3: Unemployment rate by mover type
%   Fig 4: Ergodic distributions by type (4-panel)
%   Fig 5: Unemployment rate by type + aggregate
%   Fig 6: OLF rate by type + aggregate

if ~exist('outDir','var')
    rewriteDir = fileparts(mfilename('fullpath'));
    outDir = fullfile(rewriteDir, 'Step1_FreeOmega');
end
yearsList = setdiff(1976:2024, [1977 1993]);
numYears = numel(yearsList);
Theta = 3; ThetaTot = 5;

%% NBER recession periods (start, end — calendar years, inclusive)
recessions = [1980 1980; 1981 1982; 1990 1991; 2001 2001; 2008 2009; 2020 2020];

%% Load results
omega_ts = zeros(ThetaTot, numYears);
ergMat_ts = zeros(4, Theta, numYears);
ergAgg_ts = zeros(4, numYears);

for yr = 1:numYears
    fname = fullfile(outDir, sprintf('step1_%d.mat', yr));
    if ~isfile(fname), continue; end
    tmp = load(fname); r = tmp.res;
    omega_ts(:,yr) = r.omega;
    ergMat_ts(:,:,yr) = r.ergVecs;
    ergAgg_ts(:,yr) = r.ergVecsAll;
end

%% Color scheme
c_HighN = [0.00 0.45 0.74];   % blue
c_HighU = [0.85 0.33 0.10];   % red-orange
c_HighE = [0.47 0.67 0.19];   % green
c_AllE  = [0.49 0.18 0.56];   % purple
c_AllN  = [0.64 0.08 0.18];   % dark red
c_Agg   = [0.30 0.30 0.30];   % dark gray

%% ========================================================================
%  FIGURE 1: All type shares over time
%  ========================================================================
% omega: [High-N, High-U, High-E, All-E, All-N]
fig1 = figure('Color','w','Position',[100 100 900 450]);
ax1 = axes(fig1); hold(ax1,'on');
add_recession_bars(ax1, recessions, [0 1]);
h = gobjects(5,1);
h(1) = plot(ax1, yearsList, omega_ts(4,:), '-', 'LineWidth', 2.2, 'Color', c_AllE);
h(2) = plot(ax1, yearsList, omega_ts(3,:), '-', 'LineWidth', 2.2, 'Color', c_HighE);
h(3) = plot(ax1, yearsList, omega_ts(2,:), '-', 'LineWidth', 2.2, 'Color', c_HighU);
h(4) = plot(ax1, yearsList, omega_ts(1,:), '-', 'LineWidth', 2.2, 'Color', c_HighN);
h(5) = plot(ax1, yearsList, omega_ts(5,:), '-', 'LineWidth', 2.2, 'Color', c_AllN);
legend(h, {'All-E','High-E','High-U','High-N','All-N'}, ...
    'Location','east', 'FontSize',11);
xlabel(ax1, 'Year', 'FontSize',12);
ylabel(ax1, 'Type share (\omega)', 'FontSize',12);
title(ax1, 'Step 1: Estimated type shares, free \omega', 'FontSize',14);
xlim(ax1, [1976 2024]); ylim(ax1, [0 1]);
set(ax1, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(fig1, fullfile(outDir, 'Fig1_omegas.png'), 'Resolution', 300);
savefig(fig1, fullfile(outDir, 'Fig1_omegas.fig'));

%% ========================================================================
%  FIGURE 2a-2e: Individual type shares (5 separate figures)
%  ========================================================================
typeIdx   = [5,    1,    2,    3,    4   ];
typeNames = {'All-N','High-N','High-U','High-E','All-E'};
typeCol   = {c_AllN, c_HighN, c_HighU, c_HighE, c_AllE};

for k = 1:5
    fig = figure('Color','w','Position',[100 100 500 350]);
    ax = axes(fig); hold(ax,'on');
    ydata = omega_ts(typeIdx(k),:);
    ym = max(ydata)*1.15 + 0.02;
    add_recession_bars(ax, recessions, [0 ym]);
    plot(ax, yearsList, ydata, '-', 'LineWidth', 2.2, 'Color', typeCol{k});
    title(ax, sprintf('%s share (\\omega)', typeNames{k}), 'FontSize',13);
    xlabel(ax, 'Year', 'FontSize',11);
    ylabel(ax, 'Share', 'FontSize',11);
    xlim(ax, [1976 2024]); ylim(ax, [0 ym]);
    set(ax, 'FontSize',10, 'Box','on', 'YGrid','on');
    exportgraphics(fig, fullfile(outDir, sprintf('Fig2%s_omega_%s.png', char('a'+k-1), typeNames{k})), 'Resolution', 300);
    savefig(fig, fullfile(outDir, sprintf('Fig2%s_omega_%s.fig', char('a'+k-1), typeNames{k})));
end

%% ========================================================================
%  FIGURE 3: Unemployment rate by mover type
%  ========================================================================
urate_ts = zeros(Theta, numYears);
for yr = 1:numYears
    for th = 1:Theta
        ev = ergMat_ts(:,th,yr);
        denom = ev(2)+ev(3)+ev(4);
        if denom > 0, urate_ts(th,yr) = ev(2)/denom; end
    end
end

fig3 = figure('Color','w','Position',[100 100 900 450]);
ax3 = axes(fig3); hold(ax3,'on');
add_recession_bars(ax3, recessions, [0 1]);
h3 = gobjects(3,1);
h3(1) = plot(ax3, yearsList, urate_ts(1,:), '-', 'LineWidth', 2.2, 'Color', c_HighN);
h3(2) = plot(ax3, yearsList, urate_ts(2,:), '-', 'LineWidth', 2.2, 'Color', c_HighU);
h3(3) = plot(ax3, yearsList, urate_ts(3,:), '-', 'LineWidth', 2.2, 'Color', c_HighE);
legend(h3, {'High-N','High-U','High-E'}, 'Location','northeast', 'FontSize',11);
xlabel(ax3, 'Year', 'FontSize',12);
ylabel(ax3, 'Unemployment rate', 'FontSize',12);
title(ax3, 'Step 1: Unemployment rate by mover type', 'FontSize',14);
xlim(ax3, [1976 2024]);
set(ax3, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(fig3, fullfile(outDir, 'Fig3_urate_mover.png'), 'Resolution', 300);
savefig(fig3, fullfile(outDir, 'Fig3_urate_mover.fig'));

%% ========================================================================
%  FIGURE 4: Ergodic distributions by type (4-panel)
%  ========================================================================
stateLabels = {'N (OLF)', 'U', 'E short-term', 'E long-term'};
moverNames = {'High-N', 'High-U', 'High-E'};
moverColors = {c_HighN, c_HighU, c_HighE};

fig4 = figure('Color','w','Position',[100 100 1200 800]);
for s = 1:4
    ax = subplot(2,2,s); hold(ax,'on');
    add_recession_bars(ax, recessions, [0 1]);
    h4 = gobjects(Theta,1);
    for th = 1:Theta
        h4(th) = plot(ax, yearsList, squeeze(ergMat_ts(s,th,:))', '-', ...
            'LineWidth', 2, 'Color', moverColors{th});
    end
    if s == 1, legend(h4, moverNames, 'Location','best', 'FontSize',9); end
    title(ax, sprintf('Ergodic share: %s', stateLabels{s}), 'FontSize',12);
    xlabel(ax, 'Year'); ylabel(ax, 'Share');
    xlim(ax, [1976 2024]);
    set(ax, 'FontSize',10, 'Box','on', 'YGrid','on');
end
sgtitle(fig4, 'Step 1: Ergodic distributions by type', 'FontSize',14);
exportgraphics(fig4, fullfile(outDir, 'Fig4_ergodic_by_type.png'), 'Resolution', 300);
savefig(fig4, fullfile(outDir, 'Fig4_ergodic_by_type.fig'));

%% ========================================================================
%  FIGURE 5: Unemployment rate by type + aggregate
%  ========================================================================
urate_agg = zeros(1, numYears);
for yr = 1:numYears
    ea = ergAgg_ts(:,yr);
    denom = ea(2)+ea(3)+ea(4);
    if denom > 0, urate_agg(yr) = ea(2)/denom; end
end

fig5 = figure('Color','w','Position',[100 100 900 450]);
ax5 = axes(fig5); hold(ax5,'on');
add_recession_bars(ax5, recessions, [0 1]);
h5 = gobjects(4,1);
h5(1) = plot(ax5, yearsList, urate_ts(1,:), '-', 'LineWidth', 1.8, 'Color', c_HighN);
h5(2) = plot(ax5, yearsList, urate_ts(2,:), '-', 'LineWidth', 1.8, 'Color', c_HighU);
h5(3) = plot(ax5, yearsList, urate_ts(3,:), '-', 'LineWidth', 1.8, 'Color', c_HighE);
h5(4) = plot(ax5, yearsList, urate_agg, '-', 'LineWidth', 2.5, 'Color', c_Agg);
legend(h5, {'High-N','High-U','High-E','Aggregate'}, 'Location','northeast', 'FontSize',11);
xlabel(ax5, 'Year', 'FontSize',12);
ylabel(ax5, 'Unemployment rate', 'FontSize',12);
title(ax5, 'Step 1: Unemployment rate by type and aggregate', 'FontSize',14);
xlim(ax5, [1976 2024]);
set(ax5, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(fig5, fullfile(outDir, 'Fig5_urate_all.png'), 'Resolution', 300);
savefig(fig5, fullfile(outDir, 'Fig5_urate_all.fig'));

%% ========================================================================
%  FIGURE 6: OLF rate by type + aggregate
%  ========================================================================
fig6 = figure('Color','w','Position',[100 100 900 450]);
ax6 = axes(fig6); hold(ax6,'on');
add_recession_bars(ax6, recessions, [0 1]);
h6 = gobjects(4,1);
for th = 1:Theta
    h6(th) = plot(ax6, yearsList, squeeze(ergMat_ts(1,th,:))', '-', ...
        'LineWidth', 1.8, 'Color', moverColors{th});
end
h6(4) = plot(ax6, yearsList, ergAgg_ts(1,:), '-', 'LineWidth', 2.5, 'Color', c_Agg);
legend(h6, [moverNames, {'Aggregate'}], 'Location','best', 'FontSize',11);
xlabel(ax6, 'Year', 'FontSize',12);
ylabel(ax6, 'OLF share (ergodic)', 'FontSize',12);
title(ax6, 'Step 1: OLF rate by type and aggregate', 'FontSize',14);
xlim(ax6, [1976 2024]);
set(ax6, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(fig6, fullfile(outDir, 'Fig6_olf_all.png'), 'Resolution', 300);
savefig(fig6, fullfile(outDir, 'Fig6_olf_all.fig'));

fprintf('All figures saved to %s\n', outDir);

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
