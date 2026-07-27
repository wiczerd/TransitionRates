% decomposition_table.m — Who bears the brunt of recessions?
%
% Decomposes the aggregate unemployment rate increase in each recession
% into contributions from each type.
%
% Decomposition:
%   UR_agg = sum_theta [ omega_theta * u_theta ] / LF_agg
%          = sum_theta c_theta
%   where c_theta = omega_theta * u_theta / LF_agg
%         LF_agg  = sum_theta omega_theta * (u_theta + e_theta)
%
%   Type theta's share of the brunt = Delta c_theta / Delta UR_agg

clc; clear; close all;

%% Paths
rewriteDir = fileparts(mfilename('fullpath'));
seqDir     = fullfile(rewriteDir, 'Step3_FixedOmega_Sequential');
figDir     = fullfile(seqDir, 'Figures');

%% Settings
numYears  = 47;
yearsList = setdiff(1976:2024, [1977 1993]);
ThetaTot  = 5;

typeLabels = {'High-N', 'High-U', 'High-E', 'All-E', 'All-N'};

%% Load sequential results
omega_all   = zeros(numYears, ThetaTot);
ergVecs_all = zeros(4, 4, numYears);  % 4 states x 4 types x 47 years

for yr = 1:numYears
    tmp = load(fullfile(seqDir, sprintf('step3_%d.mat', yr)));
    omega_all(yr, :)      = tmp.res.omega';
    ergVecs_all(:, :, yr) = tmp.res.ergVecs;
end

%% Compute type-level quantities for all years
% E_mass(theta, yr) = ergodic employment share (STJ+LTJ) for type theta
% U_mass(theta, yr) = ergodic unemployment share
% For All-N (type 5): U=0, E=0, N=1 by definition

E_mass = squeeze(ergVecs_all(3,:,:) + ergVecs_all(4,:,:));  % 4x47 (types 1-4)
U_mass = squeeze(ergVecs_all(2,:,:));                        % 4x47
N_mass = squeeze(ergVecs_all(1,:,:));                        % 4x47

% Add All-N as type 5: E=0, U=0
E_mass = [E_mass; zeros(1, numYears)];  % 5x47
U_mass = [U_mass; zeros(1, numYears)];
N_mass = [N_mass; ones(1, numYears)];

% Type-specific unemployment rate: U/(U+E)
Urate = U_mass ./ (U_mass + E_mass);  % 5x47
Urate(5,:) = 0;  % All-N: 0/0 -> define as 0

% Aggregate quantities for each year
LF     = zeros(1, numYears);  % aggregate labor force share
Uagg   = zeros(1, numYears);  % aggregate U mass
URagg  = zeros(1, numYears);  % aggregate U rate
c      = zeros(ThetaTot, numYears);  % type contributions to UR_agg
w      = zeros(ThetaTot, numYears);  % type labor force shares
popShare = zeros(ThetaTot, numYears); % omega (population share)

for yr = 1:numYears
    om = omega_all(yr, :);  % [High-N, High-U, High-E, All-E, All-N]

    for th = 1:ThetaTot
        Uagg(yr)  = Uagg(yr)  + om(th) * U_mass(th, yr);
        LF(yr)    = LF(yr)    + om(th) * (U_mass(th, yr) + E_mass(th, yr));
    end
    URagg(yr) = Uagg(yr) / LF(yr);

    for th = 1:ThetaTot
        c(th, yr) = om(th) * U_mass(th, yr) / LF(yr);
        w(th, yr) = om(th) * (U_mass(th, yr) + E_mass(th, yr)) / LF(yr);
        popShare(th, yr) = om(th);
    end
end

%% Define recession episodes: [peak_year, trough_year]
recessions = struct();
recessions(1).name  = '1980';
recessions(1).peak  = 1979;
recessions(1).trough = 1980;

recessions(2).name  = '1981-82';
recessions(2).peak  = 1979;
recessions(2).trough = 1982;

recessions(3).name  = '1990-91';
recessions(3).peak  = 1989;
recessions(3).trough = 1991;

recessions(4).name  = '2001';
recessions(4).peak  = 2000;
recessions(4).trough = 2001;

recessions(5).name  = '2007-09';
recessions(5).peak  = 2007;
recessions(5).trough = 2009;

recessions(6).name  = '2020';
recessions(6).peak  = 2019;
recessions(6).trough = 2020;

nRec = numel(recessions);

%% Compute decomposition for each recession
fprintf('\n');
fprintf('==========================================================================\n');
fprintf('         WHO BEARS THE BRUNT OF RECESSIONS? — Decomposition Table\n');
fprintf('==========================================================================\n');

% Store results for CSV
allRows = {};
rowIdx = 0;

for r = 1:nRec
    peakYr  = recessions(r).peak;
    troughYr = recessions(r).trough;
    ipeak   = find(yearsList == peakYr);
    itrough = find(yearsList == troughYr);

    dUR = URagg(itrough) - URagg(ipeak);

    fprintf('\n--- Recession: %s (peak=%d, trough=%d) ---\n', ...
        recessions(r).name, peakYr, troughYr);
    fprintf('Aggregate U-rate: %.2f%% -> %.2f%% (Delta = %+.2f pp)\n', ...
        URagg(ipeak)*100, URagg(itrough)*100, dUR*100);
    fprintf('\n');
    fprintf('%-8s  Pop%%   LF%%   UR_pk  UR_tr  dUR    Contrib  Share%%\n', 'Type');
    fprintf('%-8s  -----  -----  -----  -----  -----  -------  ------\n', '');

    for th = 1:ThetaTot
        dc = c(th, itrough) - c(th, ipeak);
        share_pct = dc / dUR * 100;

        fprintf('%-8s  %5.1f  %5.1f  %5.1f  %5.1f  %+5.1f  %+6.3f   %5.1f\n', ...
            typeLabels{th}, ...
            popShare(th, itrough)*100, ...
            w(th, itrough)*100, ...
            Urate(th, ipeak)*100, ...
            Urate(th, itrough)*100, ...
            (Urate(th, itrough) - Urate(th, ipeak))*100, ...
            dc*100, ...
            share_pct);

        rowIdx = rowIdx + 1;
        allRows{rowIdx, 1} = recessions(r).name;
        allRows{rowIdx, 2} = typeLabels{th};
        allRows{rowIdx, 3} = popShare(th, itrough)*100;
        allRows{rowIdx, 4} = w(th, itrough)*100;
        allRows{rowIdx, 5} = Urate(th, ipeak)*100;
        allRows{rowIdx, 6} = Urate(th, itrough)*100;
        allRows{rowIdx, 7} = (Urate(th, itrough) - Urate(th, ipeak))*100;
        allRows{rowIdx, 8} = dc*100;
        allRows{rowIdx, 9} = share_pct;
    end

    % Verify: sum of contributions = total
    dc_total = sum(c(:, itrough) - c(:, ipeak));
    fprintf('%-8s                               Total  %+6.3f   %5.1f\n', ...
        '', dc_total*100, dc_total/dUR*100);
end

%% Save CSV
T = cell2table(allRows, 'VariableNames', ...
    {'Recession','Type','PopShare_pct','LFShare_pct', ...
     'UR_peak_pct','UR_trough_pct','dUR_pp','Contribution_pp','Share_pct'});
writetable(T, fullfile(seqDir, 'decomposition_table.csv'));
fprintf('\nTable saved to: %s\n', fullfile(seqDir, 'decomposition_table.csv'));

%% =====================================================================
%  FIGURE: Stacked bar — share of brunt by type across recessions
%  =====================================================================
fig1 = figure('Position', [50 50 900 500], 'Color', 'w');

bruntShares = zeros(nRec, 4);  % exclude All-N (always 0)
recLabels = cell(nRec, 1);

for r = 1:nRec
    ipeak   = find(yearsList == recessions(r).peak);
    itrough = find(yearsList == recessions(r).trough);
    dUR = URagg(itrough) - URagg(ipeak);

    for th = 1:4
        dc = c(th, itrough) - c(th, ipeak);
        bruntShares(r, th) = dc / dUR * 100;
    end
    recLabels{r} = recessions(r).name;
end

b = bar(bruntShares, 'stacked');
b(1).FaceColor = [0.2 0.6 0.2];  % High-N
b(2).FaceColor = [0.8 0.2 0.2];  % High-U
b(3).FaceColor = [0.2 0.2 0.8];  % High-E
b(4).FaceColor = [0.5 0.5 0.5];  % All-E

set(gca, 'XTickLabel', recLabels, 'FontSize', 12);
ylabel('Share of aggregate U-rate increase (%)', 'FontSize', 12);
title('Who Bears the Brunt of Recessions?', 'FontSize', 14, 'FontWeight', 'bold');
legend({'High-N', 'High-U', 'High-E', 'All-E'}, ...
    'Location', 'eastoutside', 'FontSize', 11);
grid on; box on;
ylim([0 105]);

% Add text labels for High-U share
for r = 1:nRec
    cumul = sum(bruntShares(r, 1:1));
    ypos = cumul + bruntShares(r, 2)/2;
    text(r, ypos, sprintf('%.0f%%', bruntShares(r,2)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, ...
        'Color', 'w', 'FontWeight', 'bold');
end

saveas(fig1, fullfile(figDir, 'fig8_brunt_decomposition.png'));

%% =====================================================================
%  FIGURE: Type U-rates at peak vs trough (grouped bar)
%  =====================================================================
fig2 = figure('Position', [50 50 1000 400], 'Color', 'w');

for r = 1:nRec
    subplot(2, 3, r);
    ipeak   = find(yearsList == recessions(r).peak);
    itrough = find(yearsList == recessions(r).trough);

    barData = [Urate(1:4, ipeak)*100, Urate(1:4, itrough)*100];
    b = bar(barData);
    b(1).FaceColor = [0.6 0.8 0.6];
    b(2).FaceColor = [0.8 0.3 0.3];

    set(gca, 'XTickLabel', typeLabels(1:4), 'FontSize', 9);
    ylabel('U-rate (%)');
    title(sprintf('%s', recessions(r).name), 'FontSize', 12);
    if r == 1
        legend({sprintf('Peak (%d)', recessions(r).peak), ...
                sprintf('Trough (%d)', recessions(r).trough)}, ...
            'Location', 'northwest', 'FontSize', 8);
    else
        legend({num2str(recessions(r).peak), num2str(recessions(r).trough)}, ...
            'Location', 'northwest', 'FontSize', 8);
    end
    grid on; box on;
end

sgtitle('Type-Specific Unemployment Rates: Peak vs Trough', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(fig2, fullfile(figDir, 'fig9_urate_peak_trough.png'));

%% =====================================================================
%  COMPACT SUMMARY TABLE (for the paper)
%  =====================================================================
fprintf('\n\n');
fprintf('==========================================================================\n');
fprintf('  COMPACT TABLE: Share of Aggregate U-Rate Increase by Type (%%)\n');
fprintf('==========================================================================\n');
fprintf('%-12s  %7s  %7s  %7s  %7s  %7s\n', ...
    'Recession', 'High-N', 'High-U', 'High-E', 'All-E', 'dUR(pp)');
fprintf('%-12s  %7s  %7s  %7s  %7s  %7s\n', ...
    '----------', '------', '------', '------', '------', '------');

for r = 1:nRec
    ipeak   = find(yearsList == recessions(r).peak);
    itrough = find(yearsList == recessions(r).trough);
    dUR = URagg(itrough) - URagg(ipeak);

    shares = zeros(1, 4);
    for th = 1:4
        dc = c(th, itrough) - c(th, ipeak);
        shares(th) = dc / dUR * 100;
    end

    fprintf('%-12s  %6.1f   %6.1f   %6.1f   %6.1f   %+5.2f\n', ...
        recessions(r).name, shares(1), shares(2), shares(3), shares(4), dUR*100);
end

fprintf('\nDone.\n');
