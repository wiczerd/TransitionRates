% plot_results.m — Post-processing and plotting for Paper 2 estimation
%
% Loads per-year .mat files saved by run_estimation.m and produces figures.
% Uses named struct fields (no hardcoded row indices).

clc; clear all;

%% Load configuration
config_Paper2;

%% Which output files to load
prefix = cfg.outputPrefix;  % e.g. 'outputFO_AllE'
outDir = cfg.outputDir;

%% Type labels
if cfg.allE_estimated
    typeNames = {'High-N', 'High-U', 'High-E', 'All-E'};
else
    typeNames = {'High-N', 'High-U', 'High-E'};
end
Theta    = cfg.Theta;
ThetaTot = cfg.ThetaTot;

%% Load all years into struct arrays
numYears = cfg.numYears;
years    = cfg.yearsList;

% Pre-allocate storage
omega_ts = zeros(ThetaTot, numYears);
urate_ts = zeros(Theta, numYears);
ergMat_ts = zeros(4, Theta, numYears);
ergAgg_ts = zeros(4, numYears);
pVecs_ts  = zeros(cfg.paramsPerType, Theta, numYears);
SSR2_ts   = zeros(1, numYears);
exitflag_ts = zeros(1, numYears);

fprintf('Loading %d year files...\n', numYears);
for yr = 1:numYears
    fname = fullfile(outDir, sprintf('%s_%d.mat', prefix, yr));
    if ~isfile(fname)
        warning('Missing file: %s', fname);
        continue
    end
    tmp = load(fname);
    r = tmp.res;

    omega_ts(:, yr) = r.omega;
    urate_ts(:, yr) = r.urate(:);
    SSR2_ts(yr)     = r.SSR2;
    exitflag_ts(yr) = r.exitflag;
    pVecs_ts(:,:,yr) = r.pVecs;

    for th = 1:Theta
        ergMat_ts(:, th, yr) = r.ergVec{th};
    end
    ergAgg_ts(:, yr) = r.ergAgg;
end

sumEF1 = sum(exitflag_ts == 1);
fprintf('ExitFlag=1 in %d / %d years\n', sumEF1, numYears);

%% ========================================================================
%  Helper to build transition matrix from pVec (for computing rates)
%  ========================================================================
build_tm = @(p, isAllE) build_tm_fn(p, isAllE);

%% ========================================================================
%  FIGURE 1: Omegas over time
%  ========================================================================
f1 = figure('Color','w','Name','Omegas'); ax = axes(f1); hold(ax,'on');
allTypeNames = {'All-N', 'High-N', 'High-U', 'High-E', 'All-E'};
for k = 1:ThetaTot
    plot(ax, years, omega_ts(k,:), 'LineWidth', 2);
end
legend(ax, allTypeNames, 'Location','best');
title(ax, 'Type shares (omega) over time');
xlabel(ax, 'Year'); ylabel(ax, 'Share');
saveas(f1, fullfile(outDir, sprintf('Fig_%s_omegas.png', prefix)));

%% ========================================================================
%  FIGURE 2: Unemployment rate by type
%  ========================================================================
f2 = figure('Color','w','Name','U-rate by type'); ax = axes(f2); hold(ax,'on');
for th = 1:Theta
    plot(ax, years, urate_ts(th,:), 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'Unemployment rate by type');
xlabel(ax, 'Year'); ylabel(ax, 'U rate');
saveas(f2, fullfile(outDir, sprintf('Fig_%s_urate.png', prefix)));

%% ========================================================================
%  FIGURE 3: Aggregate ergodic distribution (N, U, ES, EL)
%  ========================================================================
f3 = figure('Color','w','Name','Aggregate ergodic'); ax = axes(f3); hold(ax,'on');
stateNames = {'N (OLF)', 'U', 'E short-term', 'E long-term'};
for s = 1:4
    plot(ax, years, ergAgg_ts(s,:), 'LineWidth', 2);
end
legend(ax, stateNames, 'Location','best');
title(ax, 'Aggregate ergodic distribution');
xlabel(ax, 'Year'); ylabel(ax, 'Share');
saveas(f3, fullfile(outDir, sprintf('Fig_%s_aggErg.png', prefix)));

%% ========================================================================
%  FIGURE 4: Ergodic distribution components by type
%  ========================================================================
ergLabels = {'OLF (state 1)', 'U (state 2)', 'E short (state 3)', 'E long (state 4)'};
for s = 1:4
    fg = figure('Color','w','Name', sprintf('Ergodic state %d', s));
    ax = axes(fg); hold(ax,'on');
    for th = 1:Theta
        plot(ax, years, squeeze(ergMat_ts(s, th, :))', 'LineWidth', 2);
    end
    legend(ax, typeNames, 'Location','best');
    title(ax, sprintf('Ergodic share: %s', ergLabels{s}));
    xlabel(ax, 'Year'); ylabel(ax, 'Share');
    saveas(fg, fullfile(outDir, sprintf('Fig_%s_erg_state%d.png', prefix, s)));
end

% Employment (states 3+4 combined)
fg = figure('Color','w','Name','Ergodic E'); ax = axes(fg); hold(ax,'on');
for th = 1:Theta
    E_share = squeeze(ergMat_ts(3,th,:) + ergMat_ts(4,th,:))';
    plot(ax, years, E_share, 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'Ergodic share: E (short + long)');
xlabel(ax, 'Year'); ylabel(ax, 'Share');
saveas(fg, fullfile(outDir, sprintf('Fig_%s_erg_E.png', prefix)));

%% ========================================================================
%  FIGURE 5: Separation rates from short-term E and long-term E
%  ========================================================================
% Sep rate from short-term E (state 3): 1 - p33 - p34 = p31 + p32
% pVec indices: p33 = idx 8, p34 = idx 9
fg_sep_ES = figure('Color','w','Name','Sep rate from ES'); ax = axes(fg_sep_ES); hold(ax,'on');
for th = 1:Theta
    sep_ES = squeeze(1 - pVecs_ts(8,th,:) - pVecs_ts(9,th,:))';
    plot(ax, years, sep_ES, 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'Separation rate from short-term E');
xlabel(ax, 'Year'); ylabel(ax, 'Rate');
saveas(fg_sep_ES, fullfile(outDir, sprintf('Fig_%s_sepRate_ES.png', prefix)));

% Sep rate from long-term E (state 4): 1 - p43 - p44 = p41 + p42
% pVec indices: p43 = idx 11, p44 = idx 12
fg_sep_EL = figure('Color','w','Name','Sep rate from EL'); ax = axes(fg_sep_EL); hold(ax,'on');
for th = 1:Theta
    sep_EL = squeeze(1 - pVecs_ts(11,th,:) - pVecs_ts(12,th,:))';
    plot(ax, years, sep_EL, 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'Separation rate from long-term E');
xlabel(ax, 'Year'); ylabel(ax, 'Rate');
saveas(fg_sep_EL, fullfile(outDir, sprintf('Fig_%s_sepRate_EL.png', prefix)));

%% ========================================================================
%  FIGURE 6: Job-finding rates
%  ========================================================================
% Job finding from U: p23 + p24 (pVec indices 5, 6)
fg_jfr_U = figure('Color','w','Name','JFR from U'); ax = axes(fg_jfr_U); hold(ax,'on');
for th = 1:Theta
    jfr = squeeze(pVecs_ts(5,th,:) + pVecs_ts(6,th,:))';
    plot(ax, years, jfr, 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'Job-finding rate from U');
xlabel(ax, 'Year'); ylabel(ax, 'Rate');
saveas(fg_jfr_U, fullfile(outDir, sprintf('Fig_%s_JFR_U.png', prefix)));

% Job finding from U to short-term E: p23 (pVec idx 5)
fg_jfr_U_ES = figure('Color','w','Name','JFR U->ES'); ax = axes(fg_jfr_U_ES); hold(ax,'on');
for th = 1:Theta
    plot(ax, years, squeeze(pVecs_ts(5,th,:))', 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'Job-finding rate from U to short-term E');
xlabel(ax, 'Year'); ylabel(ax, 'Rate');
saveas(fg_jfr_U_ES, fullfile(outDir, sprintf('Fig_%s_JFR_U_ES.png', prefix)));

% Job finding from U to long-term E: p24 (pVec idx 6)
fg_jfr_U_EL = figure('Color','w','Name','JFR U->EL'); ax = axes(fg_jfr_U_EL); hold(ax,'on');
for th = 1:Theta
    plot(ax, years, squeeze(pVecs_ts(6,th,:))', 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'Job-finding rate from U to long-term E');
xlabel(ax, 'Year'); ylabel(ax, 'Rate');
saveas(fg_jfr_U_EL, fullfile(outDir, sprintf('Fig_%s_JFR_U_EL.png', prefix)));

% Job finding from OLF (state 1): p13 + p14 (pVec indices 2, 3)
% NOTE: This fixes the bug in ProcessTimeSeries_AllE.m lines 249-250 where
% types 4 and 5 used 28*1 offset instead of 28*2 and 28*3
fg_jfr_N = figure('Color','w','Name','JFR from OLF'); ax = axes(fg_jfr_N); hold(ax,'on');
for th = 1:Theta
    jfr_N = squeeze(pVecs_ts(2,th,:) + pVecs_ts(3,th,:))';
    plot(ax, years, jfr_N, 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'Job-finding rate from OLF');
xlabel(ax, 'Year'); ylabel(ax, 'Rate');
saveas(fg_jfr_N, fullfile(outDir, sprintf('Fig_%s_JFR_N.png', prefix)));

%% ========================================================================
%  FIGURE 7: Persistence rates
%  ========================================================================
% Persistence in U: p22 (pVec idx 4)
fg_p22 = figure('Color','w','Name','Persistence U'); ax = axes(fg_p22); hold(ax,'on');
for th = 1:Theta
    plot(ax, years, squeeze(pVecs_ts(4,th,:))', 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'Persistence in U (p22)');
xlabel(ax, 'Year'); ylabel(ax, 'p22');
saveas(fg_p22, fullfile(outDir, sprintf('Fig_%s_p22.png', prefix)));

% Persistence in short-term E: p33 (pVec idx 8)
fg_p33 = figure('Color','w','Name','Persistence ES'); ax = axes(fg_p33); hold(ax,'on');
for th = 1:Theta
    plot(ax, years, squeeze(pVecs_ts(8,th,:))', 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'Persistence in short-term E (p33)');
xlabel(ax, 'Year'); ylabel(ax, 'p33');
saveas(fg_p33, fullfile(outDir, sprintf('Fig_%s_p33.png', prefix)));

% Persistence in long-term E: p44 (pVec idx 12)
fg_p44 = figure('Color','w','Name','Persistence EL'); ax = axes(fg_p44); hold(ax,'on');
for th = 1:Theta
    plot(ax, years, squeeze(pVecs_ts(12,th,:))', 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'Persistence in long-term E (p44)');
xlabel(ax, 'Year'); ylabel(ax, 'p44');
saveas(fg_p44, fullfile(outDir, sprintf('Fig_%s_p44.png', prefix)));

% Per-type persistence comparison (short vs long E)
for th = 1:Theta
    fg = figure('Color','w','Name', sprintf('Persistence %s', typeNames{th}));
    ax = axes(fg); hold(ax,'on');
    plot(ax, years, squeeze(pVecs_ts(8,th,:))', 'LineWidth', 2);
    plot(ax, years, squeeze(pVecs_ts(12,th,:))', 'LineWidth', 2);
    legend(ax, {'Short-term E (p33)', 'Long-term E (p44)'}, 'Location','best');
    title(ax, sprintf('Persistence: %s', typeNames{th}));
    xlabel(ax, 'Year'); ylabel(ax, 'Rate');
    saveas(fg, fullfile(outDir, sprintf('Fig_%s_persist_%s.png', prefix, typeNames{th})));
end

%% ========================================================================
%  FIGURE 8: Transition rates in/out of U
%  ========================================================================
% U -> N (p21 = 1 - p22 - p23 - p24)
fg_p21 = figure('Color','w','Name','U->N'); ax = axes(fg_p21); hold(ax,'on');
for th = 1:Theta
    p21 = squeeze(1 - pVecs_ts(4,th,:) - pVecs_ts(5,th,:) - pVecs_ts(6,th,:))';
    plot(ax, years, p21, 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'U to N rate (p21)');
xlabel(ax, 'Year'); ylabel(ax, 'Rate');
saveas(fg_p21, fullfile(outDir, sprintf('Fig_%s_p21.png', prefix)));

% N -> U (p12, pVec idx 1)
fg_p12 = figure('Color','w','Name','N->U'); ax = axes(fg_p12); hold(ax,'on');
for th = 1:Theta
    plot(ax, years, squeeze(pVecs_ts(1,th,:))', 'LineWidth', 2);
end
legend(ax, typeNames, 'Location','best');
title(ax, 'N to U rate (p12)');
xlabel(ax, 'Year'); ylabel(ax, 'Rate');
saveas(fg_p12, fullfile(outDir, sprintf('Fig_%s_p12.png', prefix)));

%% ========================================================================
%  Summary table
%  ========================================================================
fprintf('\n=== Summary ===\n');
fprintf('Mode: %s', cfg.mode);
if ~cfg.omegaEstimated
    fprintf(', submode: %s', cfg.fixedOmegaSubmode);
end
fprintf('\n');
fprintf('ExitFlag=1: %d / %d years\n', sumEF1, numYears);
fprintf('Figures saved to %s\n', outDir);

% =========================================================================
%  Local helper: build_tm_fn (not currently used but available for extensions)
% =========================================================================
function tm = build_tm_fn(p, isAllE)
    if isAllE
        tm = [1-p(1)-p(3), p(1), 0,    p(3); ...
              1-p(4)-p(6), p(4), 0,    p(6); ...
              0,           0,    0,    0;    ...
              1-p(10)-p(12), p(10), 0, p(12)];
    else
        tm = [1-p(1)-p(2)-p(3),  p(1), p(2), p(3); ...
              1-p(4)-p(5)-p(6),  p(4), p(5), p(6); ...
              1-p(7)-p(8)-p(9),  p(7), p(8), p(9); ...
              1-p(10)-p(11)-p(12), p(10), p(11), p(12)];
    end
end
