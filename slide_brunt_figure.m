% slide_brunt_figure.m — Clean decomposition figure for slides
% Reads existing .mat files, no re-estimation. Run time: <5 sec.
%
% Decomposition: Share_theta = omega_theta * Delta(u_theta) / sum[omega * Delta(u)]
% Sums to 100% by construction. All-N contributes zero (u=0 always).

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
ergVecs_all = zeros(4, 4, numYears);

for yr = 1:numYears
    tmp = load(fullfile(seqDir, sprintf('step3_%d.mat', yr)));
    omega_all(yr, :)      = tmp.res.omega';
    ergVecs_all(:, :, yr) = tmp.res.ergVecs;
end

%% Compute ergodic U mass by type and year
% ergVecs: 4 states x 4 types (types 1-4); state 2 = U
U_mass = squeeze(ergVecs_all(2, :, :));  % 4 x 47
U_mass = [U_mass; zeros(1, numYears)];   % 5 x 47 (All-N has u=0)

%% Define recessions
rec_def(1).name = '1980';       rec_def(1).peak = 1979; rec_def(1).trough = 1980;
rec_def(2).name = '1981-82';    rec_def(2).peak = 1979; rec_def(2).trough = 1982;
rec_def(3).name = '1990-91';    rec_def(3).peak = 1989; rec_def(3).trough = 1991;
rec_def(4).name = '2001';       rec_def(4).peak = 2000; rec_def(4).trough = 2001;
rec_def(5).name = '2007-09';    rec_def(5).peak = 2007; rec_def(5).trough = 2009;
rec_def(6).name = '2020';       rec_def(6).peak = 2019; rec_def(6).trough = 2020;
nRec = numel(rec_def);

%% Compute decomposition
% Share_theta = omega_theta * Delta(u_theta) / sum_theta'[omega_theta' * Delta(u_theta')]
shares = zeros(nRec, ThetaTot);
dUR_agg = zeros(nRec, 1);

for r = 1:nRec
    ipeak   = find(yearsList == rec_def(r).peak);
    itrough = find(yearsList == rec_def(r).trough);

    % omega is fixed across years in Phase 3, use trough year (same as peak)
    om = omega_all(itrough, :);

    % Delta u_theta for each type
    du = U_mass(:, itrough) - U_mass(:, ipeak);  % 5x1

    % Numerator: omega_theta * du_theta
    num = om(:) .* du;

    % Denominator: total
    denom = sum(num);

    shares(r, :) = (num / denom * 100)';
    dUR_agg(r) = denom * 100;  % aggregate change in U mass (pp of population)
end

%% Print summary
fprintf('\nMass decomposition: Share_theta = omega_theta * Delta(u_theta) / sum[omega*Delta(u)]\n\n');
fprintf('%-12s  %7s  %7s  %7s  %7s  %7s  %7s\n', ...
    'Recession', 'High-N', 'High-U', 'High-E', 'All-E', 'All-N', 'Sum');
for r = 1:nRec
    fprintf('%-12s  %6.1f   %6.1f   %6.1f   %6.1f   %6.1f   %6.1f\n', ...
        rec_def(r).name, shares(r,1), shares(r,2), shares(r,3), shares(r,4), shares(r,5), sum(shares(r,:)));
end

%% Figure — grouped bar, 3 mover types renormalized to 100%
recLabels = {rec_def.name};

% Renormalize: among the 3 mover types only
moverShares = shares(:, 1:3);  % High-N, High-U, High-E
for r = 1:nRec
    moverShares(r,:) = moverShares(r,:) / sum(moverShares(r,:)) * 100;
end

fprintf('\nRenormalized mover-type shares (sum to 100%%):\n');
fprintf('%-12s  %7s  %7s  %7s  %7s\n', 'Recession', 'High-N', 'High-U', 'High-E', 'Sum');
for r = 1:nRec
    fprintf('%-12s  %6.1f   %6.1f   %6.1f   %6.1f\n', ...
        rec_def(r).name, moverShares(r,1), moverShares(r,2), moverShares(r,3), sum(moverShares(r,:)));
end

fig = figure('Position', [100 100 1000 500], 'Color', 'w');

b = bar(moverShares, 'grouped');
b(1).FaceColor = [0.2 0.7 0.3];   % High-N: green
b(2).FaceColor = [0.85 0.2 0.2];  % High-U: red
b(3).FaceColor = [0.2 0.3 0.8];   % High-E: blue

set(gca, 'XTickLabel', recLabels, 'FontSize', 13, 'FontName', 'Helvetica');
ylabel('Share of mover-type \Deltau (%)', 'FontSize', 14);

legend({'High-N (~8% of pop.)', 'High-U (~4% of pop.)', 'High-E (~16% of pop.)'}, ...
    'Location', 'northwest', 'FontSize', 11);

grid on; box on;
ylim([-15 100]);
set(gca, 'YTick', -10:10:100);

% Reference line at 0
hold on;
hLine = plot(xlim, [0 0], 'k-', 'LineWidth', 0.5);
hLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold off;


%% Save
saveas(fig, fullfile(figDir, 'slide_brunt_v2.png'));
print(fig, fullfile(figDir, 'slide_brunt_v2'), '-dpng', '-r200');
fprintf('\nSaved to: %s\n', fullfile(figDir, 'slide_brunt_v2.png'));
