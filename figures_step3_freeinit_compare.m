% figures_step3_freeinit_compare.m
% Compare Phase 3 sequential (ergodic) vs free-init EM results
% Also compute fit statistics for Section 6.5

clear; clc;

%% Setup
addpath('.');
years_all = [1976, 1978:1992, 1994:2024];
Nyears = length(years_all);
recessionStarts = [1980, 1981, 1990, 2001, 2007, 2020];
recessionEnds   = [1980, 1982, 1991, 2001, 2009, 2020];

seqDir = 'Step3_FixedOmega_Sequential';
fiDir  = 'Step3_FixedOmega_FreeInit_EM';
figDir = fullfile(fiDir, 'Figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% Load data for fit stats
origDir = fullfile(fileparts(pwd), 'Replication Files and ResultsTablesFigures Files 2 (1)');
dataFile = fullfile(origDir, 'DataTimeSeriesPM_1976_2023.xlsx');
raw = readmatrix(dataFile, 'Range', 'J2:BD6562');

%% Extract results from both sets
fval_seq = zeros(Nyears,1);
fval_fi  = zeros(Nyears,1);
H_data   = zeros(Nyears,1);

% Ergodic distributions: 4 states x Nyears x 4 types (3 movers + AllE)
erg_seq = zeros(4, Nyears, 4);
erg_fi  = zeros(4, Nyears, 4);

% Init distributions (free-init only)
initdist_fi = zeros(4, Nyears, 4);

for yy = 1:Nyears
    % Sequential (ergodic)
    S = load(fullfile(seqDir, sprintf('step3_%d.mat', yy)));
    fval_seq(yy) = S.res.SSR2;

    % Free-init EM
    F = load(fullfile(fiDir, sprintf('step3_%d.mat', yy)));
    fval_fi(yy) = F.res.SSR2;

    % Ergodic distributions (already computed and saved)
    erg_seq(:, yy, :) = S.res.ergVecs;  % 4 states x 4 types
    erg_fi(:, yy, :)  = F.res.ergVecs;

    % Init dist (free-init)
    if isfield(F.res, 'initDists')
        initdist_fi(:, yy, :) = F.res.initDists;  % 4 states x 4 types
    end

    % Data entropy
    M = raw(:, yy);
    M = M / sum(M);
    idx = M > 0;
    H_data(yy) = -sum(M(idx) .* log(M(idx)));
end

%% Compute unemployment rates
urate_seq = zeros(Nyears, 4);
urate_fi  = zeros(Nyears, 4);
for yy = 1:Nyears
    for theta = 1:4
        e_s = erg_seq(3,yy,theta) + erg_seq(4,yy,theta);
        u_s = erg_seq(2,yy,theta);
        urate_seq(yy, theta) = u_s / (u_s + e_s) * 100;

        e_f = erg_fi(3,yy,theta) + erg_fi(4,yy,theta);
        u_f = erg_fi(2,yy,theta);
        urate_fi(yy, theta) = u_f / (u_f + e_f) * 100;
    end
end

%% KL divergence and fit stats
KL_seq = fval_seq - H_data;
KL_fi  = fval_fi  - H_data;
R2_seq = 1 - KL_seq ./ H_data;
R2_fi  = 1 - KL_fi  ./ H_data;

fprintf('\n=== FIT STATISTICS ===\n');
fprintf('Cross-entropy (fval):  seq mean=%.4f [%.4f, %.4f]  fi mean=%.4f [%.4f, %.4f]\n', ...
    mean(fval_seq), min(fval_seq), max(fval_seq), mean(fval_fi), min(fval_fi), max(fval_fi));
fprintf('Data entropy H(data): mean=%.4f [%.4f, %.4f]\n', mean(H_data), min(H_data), max(H_data));
fprintf('KL divergence:        seq mean=%.4f [%.4f, %.4f]  fi mean=%.4f [%.4f, %.4f]\n', ...
    mean(KL_seq), min(KL_seq), max(KL_seq), mean(KL_fi), min(KL_fi), max(KL_fi));
fprintf('R2 (1-KL/H):          seq mean=%.4f [%.4f, %.4f]  fi mean=%.4f [%.4f, %.4f]\n', ...
    mean(R2_seq), min(R2_seq), max(R2_seq), mean(R2_fi), min(R2_fi), max(R2_fi));

% Save fit stats
save(fullfile(seqDir, 'fit_stats.mat'), 'years_all', 'fval_seq', 'fval_fi', ...
    'H_data', 'KL_seq', 'KL_fi', 'R2_seq', 'R2_fi');

%% Helper for recession bands
    function add_recessions(ax, rStarts, rEnds)
        hold(ax, 'on');
        yl = ylim(ax);
        for rr = 1:length(rStarts)
            patch(ax, [rStarts(rr)-0.5, rEnds(rr)+0.5, ...
                rEnds(rr)+0.5, rStarts(rr)-0.5], ...
                [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
                'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off');
        end
        uistack(findobj(ax,'Type','patch'),'bottom');
    end

%% ---- FIGURE 1: fval improvement ----
fig1 = figure('Position',[100 100 1200 400]);
bar(years_all, fval_fi - fval_seq, 'FaceColor', [0.2 0.7 0.2]);
hold on;
yline(0, 'k-');
xlabel('Year'); ylabel('\Delta fval (free-init - sequential)');
title('Phase 3: Free Init Dist Improvement (negative = better)');
grid on;
saveas(fig1, fullfile(figDir, 'fig1_fval_improvement.png'));

%% ---- FIGURE 2: U-rate comparison (3 mover types) ----
typeNames = {'High-N', 'High-U', 'High-E'};
typeColors_seq = {[0.0 0.45 0.74], [0.85 0.33 0.1], [0.47 0.67 0.19]};
typeColors_fi  = {[0.3 0.6 0.9],   [1.0 0.6 0.4],   [0.7 0.85 0.5]};

fig2 = figure('Position',[100 100 1400 400]);
for theta = 1:3
    ax = subplot(1,3,theta);
    plot(years_all, urate_seq(:,theta), '-', 'Color', typeColors_seq{theta}, 'LineWidth', 1.5); hold on;
    plot(years_all, urate_fi(:,theta), '--', 'Color', typeColors_fi{theta}, 'LineWidth', 1.5);
    title(typeNames{theta});
    xlabel('Year'); ylabel('U-rate (%)');
    legend('Sequential (ergodic)', 'Free Init EM', 'Location', 'best');
    grid on;
    add_recessions(ax, recessionStarts, recessionEnds);
end
sgtitle('Phase 3: Unemployment Rate by Mover Type — Sequential vs Free Init');
saveas(fig2, fullfile(figDir, 'fig2_urate_comparison.png'));

%% ---- FIGURE 3: E-mass comparison ----
fig3 = figure('Position',[100 100 1400 400]);
for theta = 1:3
    ax = subplot(1,3,theta);
    emass_s = erg_seq(3,:,theta) + erg_seq(4,:,theta);
    emass_f = erg_fi(3,:,theta) + erg_fi(4,:,theta);
    plot(years_all, emass_s, '-', 'Color', typeColors_seq{theta}, 'LineWidth', 1.5); hold on;
    plot(years_all, emass_f, '--', 'Color', typeColors_fi{theta}, 'LineWidth', 1.5);
    title(typeNames{theta});
    xlabel('Year'); ylabel('Employment mass');
    legend('Sequential (ergodic)', 'Free Init EM', 'Location', 'best');
    grid on;
    add_recessions(ax, recessionStarts, recessionEnds);
end
sgtitle('Phase 3: Ergodic Employment Mass — Sequential vs Free Init');
saveas(fig3, fullfile(figDir, 'fig3_emass_comparison.png'));

%% ---- FIGURE 4: U-mass comparison ----
fig4 = figure('Position',[100 100 1400 400]);
for theta = 1:3
    ax = subplot(1,3,theta);
    plot(years_all, erg_seq(2,:,theta), '-', 'Color', typeColors_seq{theta}, 'LineWidth', 1.5); hold on;
    plot(years_all, erg_fi(2,:,theta), '--', 'Color', typeColors_fi{theta}, 'LineWidth', 1.5);
    title(typeNames{theta});
    xlabel('Year'); ylabel('U mass');
    legend('Sequential (ergodic)', 'Free Init EM', 'Location', 'best');
    grid on;
    add_recessions(ax, recessionStarts, recessionEnds);
end
sgtitle('Phase 3: Ergodic Unemployment Mass — Sequential vs Free Init');
saveas(fig4, fullfile(figDir, 'fig4_umass_comparison.png'));

%% ---- FIGURE 5: KL divergence comparison ----
fig5 = figure('Position',[100 100 1200 400]);
bar(years_all, [KL_seq, KL_fi], 'grouped');
legend('Sequential (ergodic)', 'Free Init EM', 'Location', 'northeast');
xlabel('Year'); ylabel('KL divergence');
title('KL Divergence: Model vs Data');
grid on;
saveas(fig5, fullfile(figDir, 'fig5_kl_divergence.png'));

%% ---- FIGURE 6: Init dist vs ergodic (free-init only) ----
fig6 = figure('Position',[100 100 1400 900]);
stateNames = {'OLF (s1)', 'U (s2)', 'E-short (s3)', 'E-long (s4)'};
for theta = 1:3
    for s = 1:4
        ax = subplot(3, 4, (theta-1)*4 + s);
        plot(years_all, squeeze(initdist_fi(s,:,theta)), 'b-', 'LineWidth', 1); hold on;
        plot(years_all, squeeze(erg_fi(s,:,theta)), 'r--', 'LineWidth', 1);
        if theta == 1, title(stateNames{s}); end
        if s == 1, ylabel(typeNames{theta}); end
        if theta == 3, xlabel('Year'); end
        grid on;
    end
end
legend('Init dist (free)', 'Ergodic(TM)', 'Location', 'best');
sgtitle('Phase 3 Free Init: Estimated Init Dist vs Ergodic of TM');
saveas(fig6, fullfile(figDir, 'fig6_initdist_vs_ergodic.png'));

%% ---- FIGURE 7: Noise comparison ----
noise_seq = zeros(4,1);
noise_fi  = zeros(4,1);
for theta = 1:4
    diffs_s = abs(diff(urate_seq(:,theta)));
    diffs_f = abs(diff(urate_fi(:,theta)));
    noise_seq(theta) = mean(diffs_s);
    noise_fi(theta) = mean(diffs_f);
end
fig7 = figure('Position',[100 100 600 400]);
bar(categorical({'High-N','High-U','High-E','All-E'}), [noise_seq, noise_fi], 'grouped');
legend('Sequential (ergodic)', 'Free Init EM');
ylabel('Mean |year-to-year \Delta U-rate| (pp)');
title('U-rate Noise: Sequential vs Free Init');
grid on;
saveas(fig7, fullfile(figDir, 'fig7_noise_comparison.png'));

%% Print year-by-year summary
fprintf('\n=== YEAR-BY-YEAR COMPARISON ===\n');
fprintf('Year   fval_seq   fval_fi    delta    H_data    KL_seq    KL_fi     R2_seq  R2_fi\n');
for yy = 1:Nyears
    fprintf('%d  %.4f  %.4f  %+.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
        years_all(yy), fval_seq(yy), fval_fi(yy), fval_fi(yy)-fval_seq(yy), ...
        H_data(yy), KL_seq(yy), KL_fi(yy), R2_seq(yy), R2_fi(yy));
end

fprintf('\nDone! Figures saved to %s\n', figDir);
fprintf('Fit stats saved to %s/fit_stats.mat\n', seqDir);
