% figures_step3_paper.m — Clean paper-ready figures from sequential Phase 3
%
% Generates figures for the paper draft using ONLY sequential results.
% No comparison overlays, no "Figure X:" prefixes in titles.
% Saves with paper_ prefix so existing comparison figures are preserved.

clc; clear; close all;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
step1Dir   = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
seqDir     = fullfile(rewriteDir, 'Step3_FixedOmega_Sequential');
figDir     = fullfile(seqDir, 'Figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% Settings
numYears  = 47;
yearsList = setdiff(1976:2024, [1977 1993]);
ThetaTot  = 5;

recessions = {[1980 1980], [1981 1982], [1990 1991], [2001 2001], ...
              [2008 2009], [2020 2020]};

% Colors
cHighN = [0.2 0.6 0.2];   % green
cHighU = [0.8 0.2 0.2];   % red
cHighE = [0.2 0.2 0.8];   % blue
cAllE  = [0.5 0.5 0.5];   % gray
cAllN  = [0.7 0.5 0.2];   % brown

typeLabels = {'High-N', 'High-U', 'High-E', 'All-E'};

%% Load omega_fixed (Phase 2 trends)
omf = load(fullfile(rewriteDir, 'omega_fixed.mat'));
omega_fixed = omf.omega_fixed;

%% Load Phase 1 and Sequential Phase 3
p1_omega = zeros(numYears, 5);

sq_ergVecs = zeros(4, 4, numYears);
sq_omega   = zeros(numYears, 5);

for yr = 1:numYears
    tmp1 = load(fullfile(step1Dir, sprintf('step1_%d.mat', yr)));
    p1_omega(yr, :) = tmp1.res.omega';

    tmpS = load(fullfile(seqDir, sprintf('step3_%d.mat', yr)));
    sq_ergVecs(:, :, yr) = tmpS.res.ergVecs;
    sq_omega(yr, :) = tmpS.res.omega';
end

%% Compute activity masses and unemployment rates
sq_E = squeeze(sq_ergVecs(3,:,:) + sq_ergVecs(4,:,:));  % 4x47
sq_U = squeeze(sq_ergVecs(2,:,:));
sq_N = squeeze(sq_ergVecs(1,:,:));
sq_Urate = sq_U ./ (sq_U + sq_E);

% Aggregate U-rate
sq_agg_U = zeros(1, numYears);
sq_agg_E = zeros(1, numYears);
for yr = 1:numYears
    om = sq_omega(yr,:);
    sq_agg_U(yr) = om(1)*sq_U(1,yr) + om(2)*sq_U(2,yr) + ...
                   om(3)*sq_U(3,yr) + om(4)*sq_U(4,yr);
    sq_agg_E(yr) = om(1)*sq_E(1,yr) + om(2)*sq_E(2,yr) + ...
                   om(3)*sq_E(3,yr) + om(4)*sq_E(4,yr);
end
sq_agg_Urate = sq_agg_U ./ (sq_agg_U + sq_agg_E);

%% Helper: add recession shading (HandleVisibility off)
    function add_recs(ax, recs)
        yl = ax.YLim;
        for rr = 1:numel(recs)
            rec = recs{rr};
            patch(ax, [rec(1)-0.5, rec(2)+0.5, rec(2)+0.5, rec(1)-0.5], ...
                [yl(1), yl(1), yl(2), yl(2)], ...
                [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5, ...
                'HandleVisibility','off');
        end
    end

% =====================================================================
%% PAPER FIGURE 1: Type Shares — Phase 1 estimates + Phase 2 trends
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
              [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha',0.5, ...
              'HandleVisibility','off');
    end
    plot(yearsList, p1_omega(:,i), '-o', 'Color', typeColors{i}, ...
         'MarkerSize', 3, 'MarkerFaceColor', typeColors{i}, 'LineWidth', 1);
    plot(yearsList, omega_fixed(:, trendCols(i)), '-', 'Color', 'k', ...
         'LineWidth', 2);
    title(typeNames{i}, 'FontSize', 12);
    xlim([1975 2025]);
    if i <= 3, ylabel('Population share'); end
    grid on; box on;
end
sgtitle('Type Shares: Phase 1 Estimates and Phase 2 Trends', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig1, fullfile(figDir, 'paper_fig1_omega_trends.png'));

% =====================================================================
%% PAPER FIGURE 2: Employment mass by type
% =====================================================================
fig2 = figure('Position', [50 50 800 450], 'Color', 'w');
hold on;
yl = [0 1.05]; ylim(yl);
add_recs(gca, recessions);
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
title('Employment Mass by Type', 'FontSize', 14);
xlim([1975 2025]); grid on; box on;
saveas(fig2, fullfile(figDir, 'paper_fig2_emass_bytype.png'));

% =====================================================================
%% PAPER FIGURE 3: Unemployment mass by type
% =====================================================================
fig3 = figure('Position', [50 50 800 450], 'Color', 'w');
hold on;
yl = [0 max(sq_U(:))*1.15]; ylim(yl);
add_recs(gca, recessions);
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
title('Unemployment Mass by Type', 'FontSize', 14);
xlim([1975 2025]); grid on; box on;
saveas(fig3, fullfile(figDir, 'paper_fig3_umass_bytype.png'));

% =====================================================================
%% PAPER FIGURE 4: Ergodic distributions — stacked bars, selected years
% =====================================================================
selectedYears = [1978, 1982, 1997, 2007, 2009, 2019, 2020, 2024];
selIdx = arrayfun(@(y) find(yearsList == y), selectedYears);

fig4 = figure('Position', [50 50 1000 400], 'Color', 'w');
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
sgtitle('Ergodic Distributions for Selected Years', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig4, fullfile(figDir, 'paper_fig4_ergodic_selected.png'));

% =====================================================================
%% PAPER FIGURE 5: Unemployment Rate by Type — ALL types + aggregate
% =====================================================================
fig5 = figure('Position', [50 50 1000 550], 'Color', 'w');
hold on;

yl = [0 min(100, max(sq_Urate(:))*110)];
ylim(yl);
add_recs(gca, recessions);

h1 = plot(yearsList, sq_Urate(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 2, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
h2 = plot(yearsList, sq_Urate(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 2, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
h3 = plot(yearsList, sq_Urate(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 2, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
h4 = plot(yearsList, sq_Urate(4,:)*100, '-v', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');
h5 = plot(yearsList, sq_agg_Urate*100, '-k', 'LineWidth', 2.5, 'DisplayName', 'Model aggregate');

legend([h1 h2 h3 h4 h5], 'Location', 'northeast', 'FontSize', 10);
xlabel('Year', 'FontSize', 12);
ylabel('Unemployment Rate, %', 'FontSize', 12);
title('Unemployment Rate by Type', 'FontSize', 14, 'FontWeight', 'bold');
xlim([1975 2025]); grid on; box on;
set(gca, 'FontSize', 11);
saveas(fig5, fullfile(figDir, 'paper_fig5_urate_bytype.png'));

% =====================================================================
%% PAPER FIGURE 6: U-Rate Mover Types Only
% =====================================================================
fig6 = figure('Position', [50 50 900 450], 'Color', 'w');
hold on;

yl6 = [0 max(sq_Urate(1:3,:)*100, [], 'all')*1.1];
ylim(yl6);
add_recs(gca, recessions);

plot(yearsList, sq_Urate(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 2, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, sq_Urate(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 2, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, sq_Urate(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 2, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');

legend('Location', 'northeast', 'FontSize', 10);
xlabel('Year', 'FontSize', 12);
ylabel('Unemployment Rate, %', 'FontSize', 12);
title('Unemployment Rate: Mover Types', 'FontSize', 14, 'FontWeight', 'bold');
xlim([1975 2025]); grid on; box on;
set(gca, 'FontSize', 11);
saveas(fig6, fullfile(figDir, 'paper_fig6_urate_movers.png'));

% =====================================================================
%% PAPER FIGURE 7: Brunt decomposition — stacked bar
%  (reproduced from decomposition_table.m but with paper styling)
% =====================================================================

% Load decomposition data
omega_all   = zeros(numYears, ThetaTot);
ergVecs_all = zeros(4, 4, numYears);
for yr = 1:numYears
    tmp = load(fullfile(seqDir, sprintf('step3_%d.mat', yr)));
    omega_all(yr, :) = tmp.res.omega';
    ergVecs_all(:, :, yr) = tmp.res.ergVecs;
end

E_mass = squeeze(ergVecs_all(3,:,:) + ergVecs_all(4,:,:));
U_mass = squeeze(ergVecs_all(2,:,:));
E_mass = [E_mass; zeros(1, numYears)];
U_mass = [U_mass; zeros(1, numYears)];

LF    = zeros(1, numYears);
Uagg  = zeros(1, numYears);
URagg = zeros(1, numYears);
c     = zeros(ThetaTot, numYears);

for yr = 1:numYears
    om = omega_all(yr, :);
    for th = 1:ThetaTot
        Uagg(yr) = Uagg(yr) + om(th) * U_mass(th, yr);
        LF(yr)   = LF(yr)   + om(th) * (U_mass(th, yr) + E_mass(th, yr));
    end
    URagg(yr) = Uagg(yr) / LF(yr);
    for th = 1:ThetaTot
        c(th, yr) = om(th) * U_mass(th, yr) / LF(yr);
    end
end

rec_def = struct();
rec_def(1).name = '1980';       rec_def(1).peak = 1979; rec_def(1).trough = 1980;
rec_def(2).name = '1981-82';    rec_def(2).peak = 1979; rec_def(2).trough = 1982;
rec_def(3).name = '1990-91';    rec_def(3).peak = 1989; rec_def(3).trough = 1991;
rec_def(4).name = '2001';       rec_def(4).peak = 2000; rec_def(4).trough = 2001;
rec_def(5).name = '2007-09';    rec_def(5).peak = 2007; rec_def(5).trough = 2009;
rec_def(6).name = '2020';       rec_def(6).peak = 2019; rec_def(6).trough = 2020;
nRec = numel(rec_def);

bruntShares = zeros(nRec, 4);
recLabels = cell(nRec, 1);

for r = 1:nRec
    ipeak   = find(yearsList == rec_def(r).peak);
    itrough = find(yearsList == rec_def(r).trough);
    dUR = URagg(itrough) - URagg(ipeak);
    for th = 1:4
        dc = c(th, itrough) - c(th, ipeak);
        bruntShares(r, th) = dc / dUR * 100;
    end
    recLabels{r} = rec_def(r).name;
end

fig7 = figure('Position', [50 50 900 500], 'Color', 'w');
b = bar(bruntShares, 'stacked');
b(1).FaceColor = cHighN;
b(2).FaceColor = cHighU;
b(3).FaceColor = cHighE;
b(4).FaceColor = cAllE;

set(gca, 'XTickLabel', recLabels, 'FontSize', 12);
ylabel('Share of aggregate U-rate increase (%)', 'FontSize', 12);
title('Who Bears the Brunt of Recessions?', 'FontSize', 14, 'FontWeight', 'bold');
legend({'High-N', 'High-U', 'High-E', 'All-E'}, ...
    'Location', 'eastoutside', 'FontSize', 11);
grid on; box on;
ylim([0 105]);

% Add text labels for all types
for r = 1:nRec
    cumul = 0;
    for th = 1:4
        segH = bruntShares(r, th);
        ypos = cumul + segH/2;
        % Only label if segment is large enough to fit text
        if abs(segH) > 8
            % White text on dark bars, black on light
            if th == 4  % All-E (gray bar)
                tc = 'k';
            else
                tc = 'w';
            end
            text(r, ypos, sprintf('%.0f%%', segH), ...
                'HorizontalAlignment', 'center', 'FontSize', 9, ...
                'Color', tc, 'FontWeight', 'bold');
        end
        cumul = cumul + segH;
    end
end
saveas(fig7, fullfile(figDir, 'paper_fig7_brunt_decomposition.png'));

% =====================================================================
%% PAPER FIGURE 8: Type U-rates at peak vs trough (grouped bar)
% =====================================================================
Urate_all = U_mass ./ (U_mass + E_mass);
Urate_all(5,:) = 0;
typeLabels5 = {'High-N', 'High-U', 'High-E', 'All-E', 'All-N'};

fig8 = figure('Position', [50 50 1000 400], 'Color', 'w');

for r = 1:nRec
    subplot(2, 3, r);
    ipeak   = find(yearsList == rec_def(r).peak);
    itrough = find(yearsList == rec_def(r).trough);

    barData = [Urate_all(1:4, ipeak)*100, Urate_all(1:4, itrough)*100];
    bh = bar(barData);
    bh(1).FaceColor = [0.6 0.8 0.6];
    bh(2).FaceColor = [0.8 0.3 0.3];

    set(gca, 'XTickLabel', typeLabels, 'FontSize', 9);
    ylabel('U-rate (%)');
    title(rec_def(r).name, 'FontSize', 12);
    if r == 1
        legend({sprintf('Peak (%d)', rec_def(r).peak), ...
                sprintf('Trough (%d)', rec_def(r).trough)}, ...
            'Location', 'northwest', 'FontSize', 8);
    else
        legend({num2str(rec_def(r).peak), num2str(rec_def(r).trough)}, ...
            'Location', 'northwest', 'FontSize', 8);
    end
    grid on; box on;
end

sgtitle('Type-Specific Unemployment Rates: Peak vs Trough', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig8, fullfile(figDir, 'paper_fig8_urate_peak_trough.png'));

% =====================================================================
%% LOAD TRANSITION MATRICES
% =====================================================================
Theta = 4;
tms_all = zeros(4*Theta, 4, numYears);  % 16x4x47
for yr = 1:numYears
    tmp = load(fullfile(seqDir, sprintf('step3_%d.mat', yr)));
    tms_all(:, :, yr) = tmp.res.tms;
end

% Extract key rates for each type and year
% States: 1=OLF, 2=U, 3=STJ, 4=LTJ
% TM(i,j) = probability of going from state i to state j
jobFind   = zeros(4, numYears);  % U -> E  (U->STJ + U->LTJ)
sepSTJ    = zeros(4, numYears);  % STJ -> U
sepLTJ    = zeros(4, numYears);  % LTJ -> U
UtoOLF    = zeros(4, numYears);  % U -> OLF
OLFtoU    = zeros(4, numYears);  % OLF -> U
OLFtoE    = zeros(4, numYears);  % OLF -> E  (OLF->STJ + OLF->LTJ)
EtoOLF_S  = zeros(4, numYears);  % STJ -> OLF
EtoOLF_L  = zeros(4, numYears);  % LTJ -> OLF
persOLF   = zeros(4, numYears);  % OLF -> OLF
persU     = zeros(4, numYears);  % U -> U
persSTJ   = zeros(4, numYears);  % STJ -> STJ
persLTJ   = zeros(4, numYears);  % LTJ -> LTJ
STJtoLTJ  = zeros(4, numYears);  % STJ -> LTJ (tenure)
LTJtoSTJ  = zeros(4, numYears);  % LTJ -> STJ (job switch)

for yr = 1:numYears
    for th = 1:Theta
        tm = tms_all(4*th-3:4*th, :, yr);  % 4x4
        jobFind(th,yr)  = tm(2,3) + tm(2,4);
        sepSTJ(th,yr)   = tm(3,2);
        sepLTJ(th,yr)   = tm(4,2);
        UtoOLF(th,yr)   = tm(2,1);
        OLFtoU(th,yr)   = tm(1,2);
        OLFtoE(th,yr)   = tm(1,3) + tm(1,4);
        EtoOLF_S(th,yr) = tm(3,1);
        EtoOLF_L(th,yr) = tm(4,1);
        persOLF(th,yr)  = tm(1,1);
        persU(th,yr)    = tm(2,2);
        persSTJ(th,yr)  = tm(3,3);
        persLTJ(th,yr)  = tm(4,4);
        STJtoLTJ(th,yr) = tm(3,4);
        LTJtoSTJ(th,yr) = tm(4,3);
    end
end

% =====================================================================
%% PAPER FIGURE 9: Job Finding Rate U -> E (main text)
% =====================================================================
fig9 = figure('Position', [50 50 900 450], 'Color', 'w');
hold on;
yl = [0 max(jobFind(1:3,:), [], 'all')*110]; ylim(yl);
add_recs(gca, recessions);

plot(yearsList, jobFind(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 2, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, jobFind(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 2, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, jobFind(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 2, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');

legend('Location', 'best', 'FontSize', 10);
xlabel('Year', 'FontSize', 12); ylabel('Monthly probability, %', 'FontSize', 12);
title('Job Finding Rate (U \rightarrow E)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([1975 2025]); grid on; box on;
set(gca, 'FontSize', 11);
saveas(fig9, fullfile(figDir, 'paper_fig9_job_finding.png'));

% =====================================================================
%% PAPER FIGURE 10: Job Separation Rates STJ->U and LTJ->U (main text)
% =====================================================================
fig10 = figure('Position', [50 50 1100 450], 'Color', 'w');

% Panel A: STJ -> U
subplot(1,2,1); hold on;
ylA = [0 max(sepSTJ(1:3,:), [], 'all')*120]; ylim(ylA);
add_recs(gca, recessions);
plot(yearsList, sepSTJ(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, sepSTJ(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, sepSTJ(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
plot(yearsList, sepSTJ(4,:)*100, '-v', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');
legend('Location', 'best', 'FontSize', 9);
xlabel('Year'); ylabel('Monthly probability, %');
title('(A) Short-Tenure Job Separation (STJ \rightarrow U)', 'FontSize', 12);
xlim([1975 2025]); grid on; box on;

% Panel B: LTJ -> U
subplot(1,2,2); hold on;
ylB = [0 max(sepLTJ(1:3,:), [], 'all')*120]; ylim(ylB);
add_recs(gca, recessions);
plot(yearsList, sepLTJ(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, sepLTJ(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, sepLTJ(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
plot(yearsList, sepLTJ(4,:)*100, '-v', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');
legend('Location', 'best', 'FontSize', 9);
xlabel('Year'); ylabel('Monthly probability, %');
title('(B) Long-Tenure Job Separation (LTJ \rightarrow U)', 'FontSize', 12);
xlim([1975 2025]); grid on; box on;

sgtitle('Job Separation Rates by Type', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig10, fullfile(figDir, 'paper_fig10_job_separation.png'));

% =====================================================================
%% PAPER FIGURE 11: Labor Force Transitions (main text or appendix)
%  U->OLF, OLF->U, OLF->E
% =====================================================================
fig11 = figure('Position', [50 50 1200 400], 'Color', 'w');

% Panel A: U -> OLF (discouraged workers)
subplot(1,3,1); hold on;
yl11a = [0 max(UtoOLF(1:3,:), [], 'all')*120]; ylim(yl11a);
add_recs(gca, recessions);
plot(yearsList, UtoOLF(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, UtoOLF(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, UtoOLF(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
plot(yearsList, UtoOLF(4,:)*100, '-v', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');
legend('Location', 'best', 'FontSize', 8);
xlabel('Year'); ylabel('Monthly probability, %');
title('(A) U \rightarrow OLF', 'FontSize', 11);
xlim([1975 2025]); grid on; box on;

% Panel B: OLF -> U
subplot(1,3,2); hold on;
yl11b = [0 max(OLFtoU(1:3,:), [], 'all')*120]; ylim(yl11b);
add_recs(gca, recessions);
plot(yearsList, OLFtoU(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, OLFtoU(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, OLFtoU(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
plot(yearsList, OLFtoU(4,:)*100, '-v', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');
legend('Location', 'best', 'FontSize', 8);
xlabel('Year'); ylabel('Monthly probability, %');
title('(B) OLF \rightarrow U', 'FontSize', 11);
xlim([1975 2025]); grid on; box on;

% Panel C: OLF -> E
subplot(1,3,3); hold on;
yl11c = [0 max(OLFtoE(1:3,:), [], 'all')*120]; ylim(yl11c);
add_recs(gca, recessions);
plot(yearsList, OLFtoE(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, OLFtoE(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, OLFtoE(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
plot(yearsList, OLFtoE(4,:)*100, '-v', 'Color', cAllE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cAllE, 'DisplayName', 'All-E');
legend('Location', 'best', 'FontSize', 8);
xlabel('Year'); ylabel('Monthly probability, %');
title('(C) OLF \rightarrow E', 'FontSize', 11);
xlim([1975 2025]); grid on; box on;

sgtitle('Labor Force Transition Rates by Type', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig11, fullfile(figDir, 'paper_fig11_lf_transitions.png'));

% =====================================================================
%% PAPER FIGURE 12: Persistence (diagonal elements) by type
% =====================================================================
fig12 = figure('Position', [50 50 1200 800], 'Color', 'w');

stateNames = {'OLF', 'U', 'STJ', 'LTJ'};
persData   = {persOLF, persU, persSTJ, persLTJ};

for s = 1:4
    subplot(2, 2, s); hold on;
    pd = persData{s};
    yl12 = [min(pd(1:3,:), [], 'all')*90, min(100, max(pd(1:3,:), [], 'all')*110)];
    ylim(yl12);
    add_recs(gca, recessions);

    plot(yearsList, pd(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
         'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
    plot(yearsList, pd(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 1.5, ...
         'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
    plot(yearsList, pd(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 1.5, ...
         'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');

    if s == 1, legend('Location', 'best', 'FontSize', 9); end
    xlabel('Year', 'FontSize', 10); ylabel('Monthly probability, %', 'FontSize', 10);
    title(sprintf('%s \\rightarrow %s', stateNames{s}, stateNames{s}), 'FontSize', 13);
    xlim([1975 2025]); grid on; box on;
    set(gca, 'FontSize', 10);
end
sgtitle('Persistence Rates by Type', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig12, fullfile(figDir, 'paper_fig12_persistence.png'));

% =====================================================================
%% PAPER FIGURE 13: Employment transitions — STJ<->LTJ, E->OLF
% =====================================================================
fig13 = figure('Position', [50 50 1200 800], 'Color', 'w');

% Panel A: STJ -> LTJ (tenure accumulation)
subplot(2, 2, 1); hold on;
yl13a = [0 max(STJtoLTJ(1:3,:), [], 'all')*120]; ylim(yl13a);
add_recs(gca, recessions);
plot(yearsList, STJtoLTJ(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, STJtoLTJ(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, STJtoLTJ(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
legend('Location', 'best', 'FontSize', 9);
xlabel('Year', 'FontSize', 10); ylabel('Monthly probability, %', 'FontSize', 10);
title('(A) STJ \rightarrow LTJ (tenure accumulation)', 'FontSize', 12);
xlim([1975 2025]); grid on; box on;
set(gca, 'FontSize', 10);

% Panel B: LTJ -> STJ (job switching)
subplot(2, 2, 2); hold on;
yl13b = [0 max(LTJtoSTJ(1:3,:), [], 'all')*120]; ylim(yl13b);
add_recs(gca, recessions);
plot(yearsList, LTJtoSTJ(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, LTJtoSTJ(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, LTJtoSTJ(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
legend('Location', 'best', 'FontSize', 9);
xlabel('Year', 'FontSize', 10); ylabel('Monthly probability, %', 'FontSize', 10);
title('(B) LTJ \rightarrow STJ (job switching)', 'FontSize', 12);
xlim([1975 2025]); grid on; box on;
set(gca, 'FontSize', 10);

% Panel C: STJ -> OLF (early exit)
subplot(2, 2, 3); hold on;
yl13c = [0 max(EtoOLF_S(1:3,:), [], 'all')*120]; ylim(yl13c);
add_recs(gca, recessions);
plot(yearsList, EtoOLF_S(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, EtoOLF_S(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, EtoOLF_S(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
legend('Location', 'best', 'FontSize', 9);
xlabel('Year', 'FontSize', 10); ylabel('Monthly probability, %', 'FontSize', 10);
title('(C) STJ \rightarrow OLF', 'FontSize', 12);
xlim([1975 2025]); grid on; box on;
set(gca, 'FontSize', 10);

% Panel D: LTJ -> OLF (retirement/exit)
subplot(2, 2, 4); hold on;
yl13d = [0 max(EtoOLF_L(1:3,:), [], 'all')*120]; ylim(yl13d);
add_recs(gca, recessions);
plot(yearsList, EtoOLF_L(1,:)*100, '-o', 'Color', cHighN, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighN, 'DisplayName', 'High-N');
plot(yearsList, EtoOLF_L(2,:)*100, '-s', 'Color', cHighU, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighU, 'DisplayName', 'High-U');
plot(yearsList, EtoOLF_L(3,:)*100, '-d', 'Color', cHighE, 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerFaceColor', cHighE, 'DisplayName', 'High-E');
legend('Location', 'best', 'FontSize', 9);
xlabel('Year', 'FontSize', 10); ylabel('Monthly probability, %', 'FontSize', 10);
title('(D) LTJ \rightarrow OLF', 'FontSize', 12);
xlim([1975 2025]); grid on; box on;
set(gca, 'FontSize', 10);

sgtitle('Employment Transition Rates by Type', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig13, fullfile(figDir, 'paper_fig13_employment_transitions.png'));

fprintf('Done. Paper figures saved to %s\n', figDir);
fprintf('Files: paper_fig1 through paper_fig13\n');
