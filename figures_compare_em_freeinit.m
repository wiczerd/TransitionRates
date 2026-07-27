% figures_compare_em_freeinit.m — Compare EM+FreeInit vs fmincon+FreeInit vs Ergodic
%
% Three methods:
%   1. Ergodic fmincon (best-of-both)   — blue
%   2. fmincon + Free Init Dist         — red
%   3. EM + Free Init Dist              — green
%
% Generates:
%   Fig1: Omegas comparison (5 types, 3 methods)
%   Fig2: U-rate by mover type (ergodic of TM)
%   Fig3: Fval comparison (3-method bar chart)
%   Fig4: Init dist vs ergodic — EM free-init (12 panels)
%   Fig5: Difference (initDist - ergodic) — EM free-init (12 panels)
%   Fig6: Init dist comparison: fmincon vs EM (12 panels)
%   Fig7: Fval delta EM vs fmincon free-init (head-to-head)

clc; clear;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
origDir    = fullfile(fileparts(rewriteDir), ...
               'Replication Files and ResultsTablesFigures Files 2 (1)');
bestDir    = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
fiDir      = fullfile(rewriteDir, 'Step1_FreeOmega_FreeInitDist');
emDir      = fullfile(rewriteDir, 'Step1_FreeOmega_FreeInitDist_EM');
figDir     = fullfile(emDir, 'Figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

addpath(origDir);

%% Settings
TotalYears = 47;
Theta    = 3;
ThetaTot = Theta + 2;
yearsList = setdiff(1976:2024, [1977 1993]);
map = [3, 2, 1, 1];

typeNames  = {'High-N','High-U','High-E','All-E','All-N'};
stateNames = {'OLF (s1)','U (s2)','E-short (s3)','E-long (s4)'};

%% Load all three result sets
fprintf('Loading results ...\n');

% Pre-allocate
erg_omega = nan(TotalYears, ThetaTot);  erg_fval = nan(TotalYears,1);
erg_ergVecs = nan(4, Theta, TotalYears);

fi_omega = nan(TotalYears, ThetaTot);   fi_fval = nan(TotalYears,1);
fi_ergVecs = nan(4, Theta, TotalYears);  fi_initDist = nan(4, Theta, TotalYears);

em_omega = nan(TotalYears, ThetaTot);   em_fval = nan(TotalYears,1);
em_ergVecs = nan(4, Theta, TotalYears);  em_initDist = nan(4, Theta, TotalYears);

for y = 1:TotalYears
    % Ergodic best-of-both
    tmp = load(fullfile(bestDir, sprintf('step1_%d.mat', y)));
    erg_omega(y,:)     = tmp.res.omega(:)';
    erg_ergVecs(:,:,y) = tmp.res.ergVecs;
    erg_fval(y)        = tmp.res.SSR2;

    % fmincon free-init
    tmp = load(fullfile(fiDir, sprintf('step1_%d.mat', y)));
    fi_omega(y,:)      = tmp.res.omega(:)';
    fi_ergVecs(:,:,y)  = tmp.res.ergVecs;
    fi_initDist(:,:,y) = tmp.res.initDists;
    fi_fval(y)         = tmp.res.SSR2;

    % EM free-init
    tmp = load(fullfile(emDir, sprintf('step1_%d.mat', y)));
    em_omega(y,:)      = tmp.res.omega(:)';
    em_ergVecs(:,:,y)  = tmp.res.ergVecs;
    em_initDist(:,:,y) = tmp.res.initDists;
    em_fval(y)         = tmp.res.SSR2;
end

%% Compute unemployment rates (from ergodic of TM)
erg_urate = nan(TotalYears, Theta);
fi_urate  = nan(TotalYears, Theta);
em_urate  = nan(TotalYears, Theta);
for y = 1:TotalYears
    for theta = 1:Theta
        ev = erg_ergVecs(:, theta, y);
        erg_urate(y, theta) = ev(2) / (ev(2)+ev(3)+ev(4));
        ev = fi_ergVecs(:, theta, y);
        fi_urate(y, theta)  = ev(2) / (ev(2)+ev(3)+ev(4));
        ev = em_ergVecs(:, theta, y);
        em_urate(y, theta)  = ev(2) / (ev(2)+ev(3)+ev(4));
    end
end

%% U-rates from initial distributions
fi_urate_init = nan(TotalYears, Theta);
em_urate_init = nan(TotalYears, Theta);
for y = 1:TotalYears
    for theta = 1:Theta
        pi_v = fi_initDist(:, theta, y);
        fi_urate_init(y, theta) = pi_v(2) / (pi_v(2)+pi_v(3)+pi_v(4));
        pi_v = em_initDist(:, theta, y);
        em_urate_init(y, theta) = pi_v(2) / (pi_v(2)+pi_v(3)+pi_v(4));
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

%% ===== FIG 1: OMEGA COMPARISON (3 methods) =====
fig1 = figure('Position',[100 100 1400 800]);
for k = 1:ThetaTot
    subplot(2,3,k);
    plot(yearsList, erg_omega(:,k), 'b-o', 'MarkerSize',3, 'LineWidth',1.2); hold on;
    plot(yearsList, fi_omega(:,k),  'r-s', 'MarkerSize',3, 'LineWidth',1.2);
    plot(yearsList, em_omega(:,k),  'g-^', 'MarkerSize',3, 'LineWidth',1.2);
    title(typeNames{k}, 'FontSize',12);
    xlabel('Year'); ylabel('Share');
    legend('Ergodic','fmincon FI','EM FI','Location','best');
    grid on;
    ax = gca; addRecessionShading(ax);
    ch = get(ax,'Children'); set(ax,'Children',flipud(ch));
end
sgtitle('Type Shares (Omegas): 3 Methods', 'FontSize',14);
saveas(fig1, fullfile(figDir, 'Fig1_omegas_3methods.png'));

%% ===== FIG 2: U-RATE COMPARISON (3 methods, ergodic of TM) =====
fig2 = figure('Position',[100 100 1400 400]);
for theta = 1:Theta
    subplot(1,3,theta);
    plot(yearsList, erg_urate(:,theta)*100, 'b-o', 'MarkerSize',3, 'LineWidth',1.2); hold on;
    plot(yearsList, fi_urate(:,theta)*100,  'r-s', 'MarkerSize',3, 'LineWidth',1.2);
    plot(yearsList, em_urate(:,theta)*100,  'g-^', 'MarkerSize',3, 'LineWidth',1.2);
    title(typeNames{theta}, 'FontSize',12);
    xlabel('Year'); ylabel('U rate (%)');
    legend('Ergodic','fmincon FI','EM FI','Location','best');
    grid on;
    ax = gca; addRecessionShading(ax);
    ch = get(ax,'Children'); set(ax,'Children',flipud(ch));
end
sgtitle('U Rate by Mover Type (Ergodic of TM): 3 Methods', 'FontSize',14);
saveas(fig2, fullfile(figDir, 'Fig2_urate_3methods.png'));

%% ===== FIG 3: FVAL — ALL 3 METHODS (delta from ergodic) =====
fig3 = figure('Position',[100 100 1200 500]);
delta_fi = fi_fval - erg_fval;
delta_em = em_fval - erg_fval;
bw = 0.35;
hold on;
for i = 1:TotalYears
    bar(yearsList(i) - bw/2, delta_fi(i), bw, 'FaceColor', [0.8 0.2 0.2]);
    bar(yearsList(i) + bw/2, delta_em(i), bw, 'FaceColor', [0.2 0.7 0.2]);
end
yline(0, 'k-', 'LineWidth', 1);
xlabel('Year'); ylabel('\Delta fval (method - Ergodic)');
title('Log-Likelihood Improvement vs Ergodic (negative = better)', 'FontSize',12);
legend('fmincon FI', 'EM FI', 'Location','best');
grid on;
saveas(fig3, fullfile(figDir, 'Fig3_fval_vs_ergodic.png'));

%% ===== FIG 4: INIT DIST vs ERGODIC — EM free-init (12 panels) =====
fig4 = figure('Position',[100 100 1400 900]);
for theta = 1:Theta
    for s = 1:4
        subplot(Theta, 4, (theta-1)*4 + s);
        erg_s = squeeze(em_ergVecs(s, theta, :));
        ini_s = squeeze(em_initDist(s, theta, :));
        plot(yearsList, erg_s, 'b-o', 'MarkerSize',2, 'LineWidth',1); hold on;
        plot(yearsList, ini_s, 'r-s', 'MarkerSize',2, 'LineWidth',1);
        if theta == 1, title(stateNames{s}, 'FontSize',10); end
        if s == 1, ylabel(typeNames{theta}, 'FontSize',11, 'FontWeight','bold'); end
        if theta == Theta, xlabel('Year'); end
        grid on;
        if theta == 1 && s == 1
            legend('Ergodic(tm)','Est. Init','Location','best','FontSize',7);
        end
        ax = gca; addRecessionShading(ax);
        ch = get(ax,'Children'); set(ax,'Children',flipud(ch));
    end
end
sgtitle('EM Free Init: Estimated Init Dist vs Ergodic(TM)', 'FontSize',14);
saveas(fig4, fullfile(figDir, 'Fig4_em_initdist_vs_ergodic.png'));

%% ===== FIG 5: DIFFERENCE (initDist - ergodic) — EM (12 panels) =====
fig5 = figure('Position',[100 100 1400 900]);
for theta = 1:Theta
    for s = 1:4
        subplot(Theta, 4, (theta-1)*4 + s);
        diff_s = squeeze(em_initDist(s, theta, :)) - squeeze(em_ergVecs(s, theta, :));
        bar(yearsList, diff_s, 'FaceColor', [0.4 0.6 0.8]);
        hold on; yline(0, 'k-');
        if theta == 1, title(stateNames{s}, 'FontSize',10); end
        if s == 1, ylabel(typeNames{theta}, 'FontSize',11, 'FontWeight','bold'); end
        if theta == Theta, xlabel('Year'); end
        grid on;
    end
end
sgtitle('EM Free Init: Difference (Est Init Dist - Ergodic)', 'FontSize',14);
saveas(fig5, fullfile(figDir, 'Fig5_em_initdist_minus_ergodic.png'));

%% ===== FIG 6: INIT DIST COMPARISON — fmincon vs EM (12 panels) =====
fig6 = figure('Position',[100 100 1400 900]);
for theta = 1:Theta
    for s = 1:4
        subplot(Theta, 4, (theta-1)*4 + s);
        fi_s = squeeze(fi_initDist(s, theta, :));
        em_s = squeeze(em_initDist(s, theta, :));
        plot(yearsList, fi_s, 'r-s', 'MarkerSize',2, 'LineWidth',1); hold on;
        plot(yearsList, em_s, 'g-^', 'MarkerSize',2, 'LineWidth',1);
        if theta == 1, title(stateNames{s}, 'FontSize',10); end
        if s == 1, ylabel(typeNames{theta}, 'FontSize',11, 'FontWeight','bold'); end
        if theta == Theta, xlabel('Year'); end
        grid on;
        if theta == 1 && s == 1
            legend('fmincon FI','EM FI','Location','best','FontSize',7);
        end
        ax = gca; addRecessionShading(ax);
        ch = get(ax,'Children'); set(ax,'Children',flipud(ch));
    end
end
sgtitle('Initial Distribution: fmincon vs EM (both Free Init)', 'FontSize',14);
saveas(fig6, fullfile(figDir, 'Fig6_initdist_fmincon_vs_em.png'));

%% ===== FIG 7: FVAL HEAD-TO-HEAD — EM vs fmincon free-init =====
fig7 = figure('Position',[100 100 1000 400]);
delta_em_fi = em_fval - fi_fval;
hold on;
for i = 1:TotalYears
    if delta_em_fi(i) < 0
        bar(yearsList(i), delta_em_fi(i), 'FaceColor', [0.2 0.7 0.2]);
    else
        bar(yearsList(i), delta_em_fi(i), 'FaceColor', [0.8 0.2 0.2]);
    end
end
yline(0, 'k-', 'LineWidth',1);
xlabel('Year'); ylabel('\Delta fval (EM - fmincon)');
title('EM vs fmincon (both Free Init): negative = EM wins', 'FontSize',12);
grid on;
saveas(fig7, fullfile(figDir, 'Fig7_fval_em_vs_fmincon_fi.png'));

%% ===== PRINT SUMMARY =====
fprintf('\n========================================================\n');
fprintf('SUMMARY: EM+FreeInit vs fmincon+FreeInit vs Ergodic\n');
fprintf('========================================================\n');

fprintf('\nEM vs Ergodic:    %d/%d better, %d same, %d worse\n', ...
    sum(em_fval < erg_fval - 1e-6), TotalYears, ...
    sum(abs(em_fval - erg_fval) <= 1e-6), ...
    sum(em_fval > erg_fval + 1e-6));

fprintf('fmincon FI vs Ergodic: %d/%d better, %d same, %d worse\n', ...
    sum(fi_fval < erg_fval - 1e-6), TotalYears, ...
    sum(abs(fi_fval - erg_fval) <= 1e-6), ...
    sum(fi_fval > erg_fval + 1e-6));

fprintf('EM vs fmincon FI: %d/%d better, %d same, %d worse\n', ...
    sum(em_fval < fi_fval - 1e-6), TotalYears, ...
    sum(abs(em_fval - fi_fval) <= 1e-6), ...
    sum(em_fval > fi_fval + 1e-6));

fprintf('\nMean delta (EM - Ergodic):     %.5f\n', mean(em_fval - erg_fval));
fprintf('Mean delta (fmincon - Ergodic): %.5f\n', mean(fi_fval - erg_fval));
fprintf('Mean delta (EM - fmincon FI):   %.5f\n', mean(em_fval - fi_fval));

fprintf('\nTop 10 years where EM beats fmincon FI the most:\n');
[sorted_d, si] = sort(em_fval - fi_fval);
for k = 1:min(10, TotalYears)
    fprintf('  %d: EM=%.4f, fmincon=%.4f, delta=%.4f\n', ...
        yearsList(si(k)), em_fval(si(k)), fi_fval(si(k)), sorted_d(k));
end

fprintf('\nTop 5 years where fmincon FI beats EM:\n');
for k = TotalYears:-1:max(1, TotalYears-4)
    if sorted_d(k) > 0
        fprintf('  %d: EM=%.4f, fmincon=%.4f, delta=+%.4f\n', ...
            yearsList(si(k)), em_fval(si(k)), fi_fval(si(k)), sorted_d(k));
    end
end

%% Print init dist summary stats for EM
fprintf('\n========================================================\n');
fprintf('EM Free Init: Mean abs diff (initDist - ergodic)\n');
fprintf('========================================================\n');
for theta = 1:Theta
    mad = mean(abs(squeeze(em_initDist(:,theta,:)) - squeeze(em_ergVecs(:,theta,:))), 2);
    fprintf('%s: OLF=%.4f, U=%.4f, E-short=%.4f, E-long=%.4f\n', ...
        typeNames{theta}, mad);
end

fprintf('\nCorrelation (EM initDist vs ergodic) across years:\n');
for theta = 1:Theta
    for s = 1:4
        ini_s = squeeze(em_initDist(s, theta, :));
        erg_s = squeeze(em_ergVecs(s, theta, :));
        r = corr(ini_s, erg_s);
        fprintf('  %s %s: r=%.3f\n', typeNames{theta}, stateNames{s}, r);
    end
end

fprintf('\nCorrelation (EM initDist vs fmincon initDist) across years:\n');
for theta = 1:Theta
    for s = 1:4
        em_s = squeeze(em_initDist(s, theta, :));
        fc_s = squeeze(fi_initDist(s, theta, :));
        r = corr(em_s, fc_s);
        fprintf('  %s %s: r=%.3f\n', typeNames{theta}, stateNames{s}, r);
    end
end

fprintf('\nFigures saved to %s\n', figDir);
