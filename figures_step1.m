% figures_step1.m — Publication-quality figures for Step 1 (Free Omegas)
%
% Produces:
%   Fig 1: Type shares over time (slide 51)
%   Fig 2: Individual type shares (slide 52 style, without trends)
%   Fig 3: Unemployment rate by mover type (slide 53)
%   Fig 4: Ergodic distributions by type over time
%   Fig 5: Unemployment rate by type + aggregate
%   Fig 6: OLF rate by type + aggregate

if ~exist('cfg','var')
    config_Paper2;
    cfg.mode = 'free_omega';
    cfg.Theta = 3; cfg.ThetaTot = 5;
    cfg.omegaEstimated = true; cfg.allE_estimated = false;
    cfg.outputPrefix = 'step1';
    cfg.outputDir = fullfile(cfg.dataDir, 'Step1_FreeOmega');
end

numYears = cfg.numYears;
years = cfg.yearsList;
Theta = cfg.Theta;
ThetaTot = cfg.ThetaTot;
figDir = cfg.outputDir;

%% NBER recession periods (start year, end year — inclusive, annual approx)
recessions = [1980 1980; 1981 1982; 1990 1991; 2001 2001; 2008 2009; 2020 2020];

%% Load results
omega_ts = zeros(ThetaTot, numYears);
ergMat_ts = zeros(4, Theta, numYears);
pVecs_ts = zeros(12, Theta, numYears);
ergAgg_ts = zeros(4, numYears);

for yr = 1:numYears
    fname = fullfile(cfg.outputDir, sprintf('%s_%d.mat', cfg.outputPrefix, yr));
    if ~isfile(fname), continue; end
    tmp = load(fname); r = tmp.res;
    omega_ts(:,yr) = r.omega;
    pVecs_ts(:,:,yr) = r.pVecs;
    for th = 1:Theta
        ergMat_ts(:,th,yr) = r.ergVec{th};
    end
    % Aggregate ergodic
    om = r.omega;
    ergAgg = zeros(4,1);
    for th = 1:Theta
        ergAgg = ergAgg + om(th) * r.ergVec{th};
    end
    ergAgg = ergAgg + [om(ThetaTot);0;0;om(ThetaTot-1)];
    ergAgg_ts(:,yr) = ergAgg;
end

%% Color scheme (consistent across all figures)
colors = struct();
colors.HighN  = [0.00 0.45 0.74];   % blue
colors.HighU  = [0.85 0.33 0.10];   % red-orange
colors.HighE  = [0.47 0.67 0.19];   % green
colors.AllE   = [0.49 0.18 0.56];   % purple
colors.AllN   = [0.64 0.08 0.18];   % dark red
colors.Agg    = [0.30 0.30 0.30];   % dark gray

%% ========================================================================
%  FIGURE 1: Type shares over time (Slide 51)
%  ========================================================================
% Display order: All-E, High-E, High-U, High-N, All-N (largest to smallest)
displayIdx = [4, 3, 2, 1, 5];  % indices into omega vector
dispNames  = {'All-E','High-E','High-U','High-N','All-N'};
dispColors = [colors.AllE; colors.HighE; colors.HighU; colors.HighN; colors.AllN];

fig1 = figure('Color','w','Position',[100 100 900 450]);
ax1 = axes(fig1); hold(ax1,'on');
add_recession_bars(ax1, recessions, [0 1]);
for k = 1:5
    plot(ax1, years, omega_ts(displayIdx(k),:), '-', 'LineWidth', 2.2, ...
        'Color', dispColors(k,:));
end
legend(ax1, dispNames, 'Location','east', 'FontSize',11);
xlabel(ax1, 'Year', 'FontSize',12);
ylabel(ax1, 'Type share (\omega)', 'FontSize',12);
title(ax1, 'Step 1: Estimated type shares, free \omega', 'FontSize',14);
xlim(ax1, [1976 2024]); ylim(ax1, [0 1]);
set(ax1, 'FontSize',11, 'Box','on', 'XGrid','off', 'YGrid','on');
exportgraphics(fig1, fullfile(figDir, 'Fig1_omegas.png'), 'Resolution', 300);
savefig(fig1, fullfile(figDir, 'Fig1_omegas.fig'));

%% ========================================================================
%  FIGURE 2: Individual type shares (Slide 52 style)
%  ========================================================================
fig2 = figure('Color','w','Position',[100 100 1200 350]);

% (a) High-U
ax2a = subplot(1,3,1); hold(ax2a,'on');
add_recession_bars(ax2a, recessions, [0 0.25]);
plot(ax2a, years, omega_ts(2,:), '-', 'LineWidth', 2, 'Color', colors.HighU);
title(ax2a, '(a) High-U', 'FontSize',13);
xlabel(ax2a, 'Year'); ylabel(ax2a, 'Share');
xlim(ax2a, [1976 2024]);
set(ax2a, 'FontSize',10, 'Box','on', 'YGrid','on');

% (b) High-E
ax2b = subplot(1,3,2); hold(ax2b,'on');
add_recession_bars(ax2b, recessions, [0 0.3]);
plot(ax2b, years, omega_ts(3,:), '-', 'LineWidth', 2, 'Color', colors.HighE);
title(ax2b, '(b) High-E', 'FontSize',13);
xlabel(ax2b, 'Year'); ylabel(ax2b, 'Share');
xlim(ax2b, [1976 2024]);
set(ax2b, 'FontSize',10, 'Box','on', 'YGrid','on');

% (c) All-E
ax2c = subplot(1,3,3); hold(ax2c,'on');
add_recession_bars(ax2c, recessions, [0 1]);
plot(ax2c, years, omega_ts(4,:), '-', 'LineWidth', 2, 'Color', colors.AllE);
title(ax2c, '(c) All-E', 'FontSize',13);
xlabel(ax2c, 'Year'); ylabel(ax2c, 'Share');
xlim(ax2c, [1976 2024]);
set(ax2c, 'FontSize',10, 'Box','on', 'YGrid','on');

sgtitle(fig2, 'Step 1: Type shares (individual)', 'FontSize',14);
exportgraphics(fig2, fullfile(figDir, 'Fig2_omegas_individual.png'), 'Resolution', 300);
savefig(fig2, fullfile(figDir, 'Fig2_omegas_individual.fig'));

%% ========================================================================
%  FIGURE 3: Unemployment rate by mover type (Slide 53)
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
add_recession_bars(ax3, recessions, [0 0.8]);
plot(ax3, years, urate_ts(1,:), '-', 'LineWidth', 2.2, 'Color', colors.HighN);
plot(ax3, years, urate_ts(2,:), '-', 'LineWidth', 2.2, 'Color', colors.HighU);
plot(ax3, years, urate_ts(3,:), '-', 'LineWidth', 2.2, 'Color', colors.HighE);
legend(ax3, {'High-N','High-U','High-E'}, 'Location','northeast', 'FontSize',11);
xlabel(ax3, 'Year', 'FontSize',12);
ylabel(ax3, 'Unemployment rate', 'FontSize',12);
title(ax3, 'Step 1: Unemployment rate by mover type', 'FontSize',14);
xlim(ax3, [1976 2024]);
set(ax3, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(fig3, fullfile(figDir, 'Fig3_urate_mover.png'), 'Resolution', 300);
savefig(fig3, fullfile(figDir, 'Fig3_urate_mover.fig'));

%% ========================================================================
%  FIGURE 4: Ergodic distributions by type over time
%  ========================================================================
stateLabels = {'N (OLF)', 'U', 'E short-term', 'E long-term'};
typeNames = {'High-N', 'High-U', 'High-E'};
typeColors = [colors.HighN; colors.HighU; colors.HighE];

fig4 = figure('Color','w','Position',[100 100 1200 800]);
for s = 1:4
    ax = subplot(2,2,s); hold(ax,'on');
    add_recession_bars(ax, recessions, [0 1]);
    for th = 1:Theta
        plot(ax, years, squeeze(ergMat_ts(s,th,:))', '-', ...
            'LineWidth', 2, 'Color', typeColors(th,:));
    end
    if s == 1, legend(ax, typeNames, 'Location','best', 'FontSize',9); end
    title(ax, sprintf('Ergodic share: %s', stateLabels{s}), 'FontSize',12);
    xlabel(ax, 'Year'); ylabel(ax, 'Share');
    xlim(ax, [1976 2024]); ylim(ax, [0 max(squeeze(max(max(ergMat_ts(s,:,:)))))*1.1+0.01]);
    set(ax, 'FontSize',10, 'Box','on', 'YGrid','on');
end
sgtitle(fig4, 'Step 1: Ergodic distributions by type', 'FontSize',14);
exportgraphics(fig4, fullfile(figDir, 'Fig4_ergodic_by_type.png'), 'Resolution', 300);
savefig(fig4, fullfile(figDir, 'Fig4_ergodic_by_type.fig'));

%% ========================================================================
%  FIGURE 5: Unemployment rate by type + aggregate
%  ========================================================================
% Aggregate U-rate
urate_agg = zeros(1, numYears);
for yr = 1:numYears
    ea = ergAgg_ts(:,yr);
    denom = ea(2)+ea(3)+ea(4);
    if denom > 0, urate_agg(yr) = ea(2)/denom; end
end

fig5 = figure('Color','w','Position',[100 100 900 450]);
ax5 = axes(fig5); hold(ax5,'on');
add_recession_bars(ax5, recessions, [0 0.8]);
plot(ax5, years, urate_ts(1,:), '-', 'LineWidth', 1.8, 'Color', colors.HighN);
plot(ax5, years, urate_ts(2,:), '-', 'LineWidth', 1.8, 'Color', colors.HighU);
plot(ax5, years, urate_ts(3,:), '-', 'LineWidth', 1.8, 'Color', colors.HighE);
plot(ax5, years, urate_agg, '-', 'LineWidth', 2.5, 'Color', colors.Agg);
legend(ax5, {'High-N','High-U','High-E','Aggregate'}, 'Location','northeast', 'FontSize',11);
xlabel(ax5, 'Year', 'FontSize',12);
ylabel(ax5, 'Unemployment rate', 'FontSize',12);
title(ax5, 'Step 1: Unemployment rate by type and aggregate', 'FontSize',14);
xlim(ax5, [1976 2024]);
set(ax5, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(fig5, fullfile(figDir, 'Fig5_urate_all.png'), 'Resolution', 300);
savefig(fig5, fullfile(figDir, 'Fig5_urate_all.fig'));

%% ========================================================================
%  FIGURE 6: OLF rate by type + aggregate
%  ========================================================================
% OLF rate = ergodic mass in state 1 (already a share of population)
olf_agg = ergAgg_ts(1,:);

fig6 = figure('Color','w','Position',[100 100 900 450]);
ax6 = axes(fig6); hold(ax6,'on');
add_recession_bars(ax6, recessions, [0 1]);
for th = 1:Theta
    plot(ax6, years, squeeze(ergMat_ts(1,th,:))', '-', ...
        'LineWidth', 1.8, 'Color', typeColors(th,:));
end
plot(ax6, years, olf_agg, '-', 'LineWidth', 2.5, 'Color', colors.Agg);
legend(ax6, [typeNames, {'Aggregate'}], 'Location','best', 'FontSize',11);
xlabel(ax6, 'Year', 'FontSize',12);
ylabel(ax6, 'OLF share (ergodic)', 'FontSize',12);
title(ax6, 'Step 1: OLF rate by type and aggregate', 'FontSize',14);
xlim(ax6, [1976 2024]);
set(ax6, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(fig6, fullfile(figDir, 'Fig6_olf_all.png'), 'Resolution', 300);
savefig(fig6, fullfile(figDir, 'Fig6_olf_all.fig'));

fprintf('All figures saved to %s\n', figDir);

%% ========================================================================
%  Helper: add recession shading
%  ========================================================================
function add_recession_bars(ax, recessions, ylims)
    for r = 1:size(recessions,1)
        x1 = recessions(r,1) - 0.5;
        x2 = recessions(r,2) + 0.5;
        fill(ax, [x1 x2 x2 x1], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
            [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.5);
    end
    % Send recession bars to back
    children = get(ax, 'Children');
    nRec = size(recessions,1);
    if numel(children) >= nRec
        set(ax, 'Children', [children(nRec+1:end); children(1:nRec)]);
    end
end
