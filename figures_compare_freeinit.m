% figures_compare_freeinit.m — Compare free-init-dist vs ergodic results
%
% Generates:
%   Fig1: Omegas comparison (5 types, 2 methods)
%   Fig2: U-rate by mover type comparison
%   Fig3: Fval comparison (bar chart)
%   Fig4: Estimated initial dist vs ergodic dist for each mover type
%   Fig5: Difference (initDist - ergodic) over time, by type and state
%   Table: Initial distributions for selected years

clc; clear;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(fileparts(rewriteDir), ...
               'Replication Files and ResultsTablesFigures Files 2 (1)');
bestDir    = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
fiDir      = fullfile(rewriteDir, 'Step1_FreeOmega_FreeInitDist');
figDir     = fullfile(fiDir, 'Figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

addpath(origDir);

%% Settings
TotalYears = 47;
Theta    = 3;
ThetaTot = Theta + 2;
yearsList = setdiff(1976:2024, [1977 1993]);
map = [3, 2, 1, 1];  % state -> activity: N, U, E, E

typeNames  = {'High-N','High-U','High-E','All-E','All-N'};
stateNames = {'OLF (s1)','U (s2)','E-short (s3)','E-long (s4)'};
recYears   = [1980 1981 1982 1990 1991 2001 2008 2009 2010 2020];

%% Load both result sets
fprintf('Loading results ...\n');
fi_omega    = nan(TotalYears, ThetaTot);
fi_ergVecs  = nan(4, Theta, TotalYears);
fi_initDist = nan(4, Theta, TotalYears);
fi_fval     = nan(TotalYears, 1);
fi_pVecs    = nan(12, Theta, TotalYears);

erg_omega   = nan(TotalYears, ThetaTot);
erg_ergVecs = nan(4, Theta, TotalYears);
erg_fval    = nan(TotalYears, 1);
erg_pVecs   = nan(12, Theta, TotalYears);

for y = 1:TotalYears
    % Free init dist
    tmp = load(fullfile(fiDir, sprintf('step1_%d.mat', y)));
    fi_omega(y,:)       = tmp.res.omega(:)';
    fi_ergVecs(:,:,y)   = tmp.res.ergVecs;
    fi_initDist(:,:,y)  = tmp.res.initDists;
    fi_fval(y)          = tmp.res.SSR2;
    fi_pVecs(:,:,y)     = tmp.res.pVecs;

    % Ergodic best-of-both
    tmp2 = load(fullfile(bestDir, sprintf('step1_%d.mat', y)));
    erg_omega(y,:)      = tmp2.res.omega(:)';
    erg_ergVecs(:,:,y)  = tmp2.res.ergVecs;
    erg_fval(y)         = tmp2.res.SSR2;
    erg_pVecs(:,:,y)    = tmp2.res.pVecs;
end

%% Compute unemployment rates by type (from ergodic distributions)
% U rate = ergVec(2) / (ergVec(2) + ergVec(3) + ergVec(4))  [excluding OLF]
% Or simpler: U rate = state2 / (1 - state1)  among labor force participants
fi_urate  = nan(TotalYears, Theta);
erg_urate = nan(TotalYears, Theta);
for y = 1:TotalYears
    for theta = 1:Theta
        ev_fi  = fi_ergVecs(:, theta, y);
        ev_erg = erg_ergVecs(:, theta, y);
        fi_urate(y, theta)  = ev_fi(2) / (ev_fi(2) + ev_fi(3) + ev_fi(4));
        erg_urate(y, theta) = ev_erg(2) / (ev_erg(2) + ev_erg(3) + ev_erg(4));
    end
end

%% Compute initial-dist-implied unemployment rates
fi_urate_init = nan(TotalYears, Theta);
for y = 1:TotalYears
    for theta = 1:Theta
        pi_v = fi_initDist(:, theta, y);
        fi_urate_init(y, theta) = pi_v(2) / (pi_v(2) + pi_v(3) + pi_v(4));
    end
end

%% Recession shading helper
recessionPeriods = [1980 1982; 1990 1991; 2001 2001; 2008 2009; 2020 2020];
addRecessionShading = @(ax) arrayfun(@(i) ...
    patch(ax, [recessionPeriods(i,1)-0.5 recessionPeriods(i,2)+0.5 ...
               recessionPeriods(i,2)+0.5 recessionPeriods(i,1)-0.5], ...
          [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5), ...
    1:size(recessionPeriods,1));

%% ===== FIG 1: OMEGA COMPARISON =====
fig1 = figure('Position',[100 100 1400 800]);
colors = lines(5);
for k = 1:ThetaTot
    subplot(2,3,k);
    plot(yearsList, erg_omega(:,k), 'b-o', 'MarkerSize',3, 'LineWidth',1.2); hold on;
    plot(yearsList, fi_omega(:,k), 'r-s', 'MarkerSize',3, 'LineWidth',1.2);
    title(typeNames{k}, 'FontSize',12);
    xlabel('Year'); ylabel('Share');
    legend('Ergodic','Free Init','Location','best');
    grid on;
    ax = gca;
    addRecessionShading(ax);
    % Bring lines to front
    ch = get(ax,'Children');
    set(ax,'Children',flipud(ch));
end
sgtitle('Type Shares (Omegas): Ergodic vs Free Init Dist', 'FontSize',14);
saveas(fig1, fullfile(figDir, 'Fig1_omegas_comparison.png'));

%% ===== FIG 2: U-RATE COMPARISON =====
fig2 = figure('Position',[100 100 1400 400]);
for theta = 1:Theta
    subplot(1,3,theta);
    plot(yearsList, erg_urate(:,theta)*100, 'b-o', 'MarkerSize',3, 'LineWidth',1.2); hold on;
    plot(yearsList, fi_urate(:,theta)*100, 'r-s', 'MarkerSize',3, 'LineWidth',1.2);
    title(typeNames{theta}, 'FontSize',12);
    xlabel('Year'); ylabel('U rate (%)');
    legend('Ergodic','Free Init','Location','best');
    grid on;
    ax = gca;
    addRecessionShading(ax);
    ch = get(ax,'Children');
    set(ax,'Children',flipud(ch));
end
sgtitle('Unemployment Rate by Mover Type (from Ergodic Dist of TM)', 'FontSize',14);
saveas(fig2, fullfile(figDir, 'Fig2_urate_comparison.png'));

%% ===== FIG 3: FVAL COMPARISON =====
fig3 = figure('Position',[100 100 1000 400]);
delta = fi_fval - erg_fval;
bar(yearsList, delta);
hold on;
yline(0, 'k-', 'LineWidth',1);
xlabel('Year'); ylabel('\Delta fval (Free Init - Ergodic)');
title('Log-Likelihood Improvement: Free Init vs Ergodic (negative = better)', 'FontSize',12);
grid on;
% Color bars: green for improvement, red for worse
for i = 1:length(delta)
    if delta(i) < 0
        bar(yearsList(i), delta(i), 'FaceColor', [0.2 0.7 0.2]);
    else
        bar(yearsList(i), delta(i), 'FaceColor', [0.8 0.2 0.2]);
    end
    hold on;
end
yline(0, 'k-', 'LineWidth',1);
saveas(fig3, fullfile(figDir, 'Fig3_fval_delta.png'));

%% ===== FIG 4: INITIAL DIST vs ERGODIC for each mover type =====
% For each type, plot 4 states over time: estimated initDist vs ergodic(tm)
fig4 = figure('Position',[100 100 1400 900]);
for theta = 1:Theta
    for s = 1:4
        subplot(Theta, 4, (theta-1)*4 + s);
        erg_s = squeeze(fi_ergVecs(s, theta, :));  % ergodic from free-init tm
        ini_s = squeeze(fi_initDist(s, theta, :)); % estimated init dist
        plot(yearsList, erg_s, 'b-o', 'MarkerSize',2, 'LineWidth',1); hold on;
        plot(yearsList, ini_s, 'r-s', 'MarkerSize',2, 'LineWidth',1);
        if theta == 1
            title(stateNames{s}, 'FontSize',10);
        end
        if s == 1
            ylabel(typeNames{theta}, 'FontSize',11, 'FontWeight','bold');
        end
        if theta == Theta
            xlabel('Year');
        end
        grid on;
        if theta == 1 && s == 1
            legend('Ergodic(tm)','Est. Init','Location','best','FontSize',7);
        end
        ax = gca;
        addRecessionShading(ax);
        ch = get(ax,'Children');
        set(ax,'Children',flipud(ch));
    end
end
sgtitle('Estimated Initial Distribution vs Ergodic(TM) by Type and State', 'FontSize',14);
saveas(fig4, fullfile(figDir, 'Fig4_initdist_vs_ergodic.png'));

%% ===== FIG 5: DIFFERENCE (initDist - ergodic) =====
fig5 = figure('Position',[100 100 1400 900]);
for theta = 1:Theta
    for s = 1:4
        subplot(Theta, 4, (theta-1)*4 + s);
        diff_s = squeeze(fi_initDist(s, theta, :)) - squeeze(fi_ergVecs(s, theta, :));
        bar(yearsList, diff_s, 'FaceColor', [0.4 0.6 0.8]);
        hold on; yline(0, 'k-');
        if theta == 1
            title(stateNames{s}, 'FontSize',10);
        end
        if s == 1
            ylabel(typeNames{theta}, 'FontSize',11, 'FontWeight','bold');
        end
        if theta == Theta
            xlabel('Year');
        end
        grid on;
    end
end
sgtitle('Difference: Estimated Init Dist - Ergodic(TM)', 'FontSize',14);
saveas(fig5, fullfile(figDir, 'Fig5_initdist_minus_ergodic.png'));

%% ===== FIG 6: U-rate from init dist vs ergodic dist =====
fig6 = figure('Position',[100 100 1400 400]);
for theta = 1:Theta
    subplot(1,3,theta);
    plot(yearsList, erg_urate(:,theta)*100, 'b-o', 'MarkerSize',3, 'LineWidth',1.2); hold on;
    plot(yearsList, fi_urate(:,theta)*100, 'r-s', 'MarkerSize',3, 'LineWidth',1.2);
    plot(yearsList, fi_urate_init(:,theta)*100, 'g-^', 'MarkerSize',3, 'LineWidth',1.2);
    title(typeNames{theta}, 'FontSize',12);
    xlabel('Year'); ylabel('U rate (%)');
    legend('Ergodic run','Free-init (erg of tm)','Free-init (init dist)','Location','best');
    grid on;
    ax = gca;
    addRecessionShading(ax);
    ch = get(ax,'Children');
    set(ax,'Children',flipud(ch));
end
sgtitle('U Rate: Ergodic(TM) vs Estimated Initial Distribution', 'FontSize',14);
saveas(fig6, fullfile(figDir, 'Fig6_urate_initdist.png'));

%% ===== TABLE: Initial distributions for selected years =====
fprintf('\n========================================================\n');
fprintf('INITIAL DISTRIBUTIONS vs ERGODIC for SELECTED YEARS\n');
fprintf('========================================================\n');
selectedYears = [1978, 1982, 2001, 2007, 2008, 2009, 2010, 2019, 2020, 2021];

for yi = 1:length(selectedYears)
    yr = selectedYears(yi);
    idx = find(yearsList == yr);
    if isempty(idx), continue; end

    isRec = any(yr == recYears);
    recTag = '';
    if isRec, recTag = ' [RECESSION]'; end

    fprintf('\n--- Year %d%s (fval: erg=%.4f, free=%.4f, delta=%.4f) ---\n', ...
        yr, recTag, erg_fval(idx), fi_fval(idx), fi_fval(idx)-erg_fval(idx));
    fprintf('  Omegas (erg): '); fprintf('%.3f ', erg_omega(idx,:)); fprintf('\n');
    fprintf('  Omegas (fi):  '); fprintf('%.3f ', fi_omega(idx,:)); fprintf('\n');

    for theta = 1:Theta
        fprintf('  %s:\n', typeNames{theta});
        erg_v = fi_ergVecs(:, theta, idx);
        ini_v = fi_initDist(:, theta, idx);
        fprintf('    Ergodic(tm): [%.4f  %.4f  %.4f  %.4f]  U-rate=%.1f%%\n', ...
            erg_v, erg_v(2)/(erg_v(2)+erg_v(3)+erg_v(4))*100);
        fprintf('    Est. Init:   [%.4f  %.4f  %.4f  %.4f]  U-rate=%.1f%%\n', ...
            ini_v, ini_v(2)/(ini_v(2)+ini_v(3)+ini_v(4))*100);
        fprintf('    Diff:        [%+.4f %+.4f %+.4f %+.4f]\n', ini_v - erg_v);
    end
end

%% Summary statistics
fprintf('\n========================================================\n');
fprintf('SUMMARY: Mean absolute difference (initDist - ergodic)\n');
fprintf('========================================================\n');
for theta = 1:Theta
    mad = mean(abs(squeeze(fi_initDist(:,theta,:)) - squeeze(fi_ergVecs(:,theta,:))), 2);
    fprintf('%s: OLF=%.4f, U=%.4f, E-short=%.4f, E-long=%.4f\n', ...
        typeNames{theta}, mad);
end

% Correlation between init dist and ergodic
fprintf('\n--- Correlation (initDist vs ergodic) across years ---\n');
for theta = 1:Theta
    for s = 1:4
        ini_series = squeeze(fi_initDist(s, theta, :));
        erg_series = squeeze(fi_ergVecs(s, theta, :));
        r = corr(ini_series, erg_series);
        fprintf('%s %s: r=%.3f\n', typeNames{theta}, stateNames{s}, r);
    end
end

fprintf('\nFigures saved to %s\n', figDir);
