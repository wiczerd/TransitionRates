% figures_compare_em.m — Compare EM vs fmincon best-of-both results
%
% Produces:
%   FigE1: Type shares — fmincon (dashed) vs EM (solid), 5 panels
%   FigE2: U-rate by mover type — fmincon vs EM, 3 panels
%   FigE3: Fval comparison bar chart
%   FigE4: EM results — all type shares
%   FigE5: EM results — u-rate by mover type

rewriteDir = fileparts(mfilename('fullpath'));
bestDir = fullfile(rewriteDir, 'Step1_FreeOmega_Best');
emDir   = fullfile(rewriteDir, 'Step1_FreeOmega_EM');
figDir  = fullfile(emDir, 'Figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

yearsList = setdiff(1976:2024, [1977 1993]);
numYears  = numel(yearsList);
Theta = 3; ThetaTot = 5;
recessions = [1980 1980; 1981 1982; 1990 1991; 2001 2001; 2008 2009; 2020 2020];

%% Load both sets
omega_fc  = zeros(ThetaTot, numYears);
omega_em  = zeros(ThetaTot, numYears);
erg_fc    = zeros(4, Theta, numYears);
erg_em    = zeros(4, Theta, numYears);
ergAgg_fc = zeros(4, numYears);
ergAgg_em = zeros(4, numYears);
fval_fc   = nan(1, numYears);
fval_em   = nan(1, numYears);

for yr = 1:numYears
    tmp = load(fullfile(bestDir, sprintf('step1_%d.mat', yr)));
    r = tmp.res;
    omega_fc(:,yr) = r.omega; erg_fc(:,:,yr) = r.ergVecs;
    ergAgg_fc(:,yr) = r.ergVecsAll; fval_fc(yr) = r.SSR2;

    tmp = load(fullfile(emDir, sprintf('step1_%d.mat', yr)));
    r = tmp.res;
    omega_em(:,yr) = r.omega; erg_em(:,:,yr) = r.ergVecs;
    ergAgg_em(:,yr) = r.ergVecsAll; fval_em(yr) = r.SSR2;
end

%% Colors
c_HighN = [0.00 0.45 0.74];
c_HighU = [0.85 0.33 0.10];
c_HighE = [0.47 0.67 0.19];
c_AllE  = [0.49 0.18 0.56];
c_AllN  = [0.64 0.08 0.18];
typeIdx   = [4, 3, 2, 1, 5];
typeNames = {'All-E','High-E','High-U','High-N','All-N'};
typeCol   = {c_AllE, c_HighE, c_HighU, c_HighN, c_AllN};
moverNames = {'High-N','High-U','High-E'};
moverCol   = {c_HighN, c_HighU, c_HighE};

%% U-rates
urate_fc = zeros(Theta, numYears);
urate_em = zeros(Theta, numYears);
for yr = 1:numYears
    for th = 1:Theta
        ev = erg_fc(:,th,yr); d = ev(2)+ev(3)+ev(4);
        if d > 0, urate_fc(th,yr) = ev(2)/d; end
        ev = erg_em(:,th,yr); d = ev(2)+ev(3)+ev(4);
        if d > 0, urate_em(th,yr) = ev(2)/d; end
    end
end

%% FigE1: Type shares comparison — 5 panels
figE1 = figure('Color','w','Position',[50 50 1400 900]);
for k = 1:5
    ax = subplot(3,2,k); hold(ax,'on');
    idx = typeIdx(k);
    ym = max(max(omega_fc(idx,:)), max(omega_em(idx,:)))*1.15 + 0.02;
    add_recession_bars(ax, recessions, [0 ym]);
    plot(ax, yearsList, omega_fc(idx,:), '--', 'LineWidth', 1.5, 'Color', [0.6 0.6 0.6]);
    plot(ax, yearsList, omega_em(idx,:), '-',  'LineWidth', 2.0, 'Color', typeCol{k});
    title(ax, typeNames{k}, 'FontSize', 12);
    xlim(ax, [1976 2024]); ylim(ax, [0 ym]);
    set(ax, 'FontSize', 10, 'Box', 'on', 'YGrid', 'on');
    if k == 1, legend(ax, {'fmincon','EM'}, 'Location','best', 'FontSize', 9); end
    if k >= 4, xlabel(ax, 'Year'); end
    ylabel(ax, '\omega');
end
ax6 = subplot(3,2,6); axis(ax6, 'off');
nEMbetter = sum(fval_em < fval_fc - 1e-6);
nFCbetter = sum(fval_fc < fval_em - 1e-6);
text(ax6, 0.05, 0.8, sprintf('EM better: %d years', nEMbetter), 'FontSize', 12, 'Color', [0.2 0.6 0.2]);
text(ax6, 0.05, 0.6, sprintf('fmincon better: %d years', nFCbetter), 'FontSize', 12, 'Color', [0.8 0.2 0.2]);
text(ax6, 0.05, 0.4, sprintf('Agree (|delta|<1e-4): %d years', sum(abs(fval_fc-fval_em)<1e-4)), 'FontSize', 12);
sgtitle(figE1, 'Type Shares: fmincon best-of-both (dashed) vs EM+refine (solid)', 'FontSize', 14);
exportgraphics(figE1, fullfile(figDir, 'FigE1_omegas_em_vs_fc.png'), 'Resolution', 300);

%% FigE2: U-rate comparison — 3 panels
figE2 = figure('Color','w','Position',[100 100 1200 500]);
for th = 1:Theta
    ax = subplot(1,3,th); hold(ax,'on');
    add_recession_bars(ax, recessions, [0 1]);
    plot(ax, yearsList, urate_fc(th,:), '--', 'LineWidth', 1.5, 'Color', [0.6 0.6 0.6]);
    plot(ax, yearsList, urate_em(th,:), '-',  'LineWidth', 2.0, 'Color', moverCol{th});
    title(ax, moverNames{th}, 'FontSize', 12);
    xlim(ax, [1976 2024]); ylim(ax, [0 1]);
    set(ax, 'FontSize', 10, 'Box', 'on', 'YGrid', 'on');
    xlabel(ax, 'Year'); ylabel(ax, 'U-rate');
    if th == 1, legend(ax, {'fmincon','EM'}, 'Location','best', 'FontSize', 9); end
end
sgtitle(figE2, 'Unemployment Rate: fmincon vs EM', 'FontSize', 14);
exportgraphics(figE2, fullfile(figDir, 'FigE2_urate_em_vs_fc.png'), 'Resolution', 300);

%% FigE3: Fval comparison
figE3 = figure('Color','w','Position',[100 100 900 400]);
ax = axes(figE3); hold(ax,'on');
delta = fval_fc - fval_em;  % positive = EM better
bar_colors = zeros(numYears, 3);
for i = 1:numYears
    if delta(i) > 1e-6
        bar_colors(i,:) = [0.2 0.6 0.2];  % green = EM better
    else
        bar_colors(i,:) = [0.8 0.2 0.2];  % red = fmincon better
    end
end
for i = 1:numYears
    bar(ax, yearsList(i), delta(i), 0.7, 'FaceColor', bar_colors(i,:), 'EdgeColor', 'none');
end
yline(ax, 0, 'k-', 'LineWidth', 0.5);
xlabel(ax, 'Year', 'FontSize', 12);
ylabel(ax, '\Delta fval (fmincon - EM)', 'FontSize', 12);
title(ax, 'EM vs fmincon: positive = EM found better solution', 'FontSize', 13);
set(ax, 'FontSize', 10, 'Box', 'on', 'YGrid', 'on');
xlim(ax, [1975 2025]);
exportgraphics(figE3, fullfile(figDir, 'FigE3_fval_em_vs_fc.png'), 'Resolution', 300);

%% FigE4: EM standalone — type shares
figE4 = figure('Color','w','Position',[100 100 900 450]);
ax = axes(figE4); hold(ax,'on');
add_recession_bars(ax, recessions, [0 1]);
h = gobjects(5,1);
for k = 1:5
    h(k) = plot(ax, yearsList, omega_em(typeIdx(k),:), '-', 'LineWidth', 2.2, 'Color', typeCol{k});
end
legend(h, typeNames, 'Location','east', 'FontSize',11);
xlabel(ax, 'Year', 'FontSize',12);
ylabel(ax, 'Type share (\omega)', 'FontSize',12);
title(ax, 'Step 1: Estimated type shares (EM algorithm)', 'FontSize',14);
xlim(ax, [1976 2024]); ylim(ax, [0 1]);
set(ax, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(figE4, fullfile(figDir, 'FigE4_omegas_em.png'), 'Resolution', 300);

%% FigE5: EM standalone — u-rate
figE5 = figure('Color','w','Position',[100 100 900 450]);
ax = axes(figE5); hold(ax,'on');
add_recession_bars(ax, recessions, [0 1]);
h5 = gobjects(3,1);
for th = 1:Theta
    h5(th) = plot(ax, yearsList, urate_em(th,:), '-', 'LineWidth', 2.2, 'Color', moverCol{th});
end
legend(h5, moverNames, 'Location','northeast', 'FontSize',11);
xlabel(ax, 'Year', 'FontSize',12);
ylabel(ax, 'Unemployment rate', 'FontSize',12);
title(ax, 'Step 1: Unemployment rate by mover type (EM algorithm)', 'FontSize',14);
xlim(ax, [1976 2024]);
set(ax, 'FontSize',11, 'Box','on', 'YGrid','on');
exportgraphics(figE5, fullfile(figDir, 'FigE5_urate_em.png'), 'Resolution', 300);

fprintf('\nAll figures saved to:\n  %s\n', figDir);
fprintf('\nAgreement summary:\n');
fprintf('  Years where |delta fval| < 0.001: %d / %d\n', sum(abs(fval_fc - fval_em) < 0.001), numYears);
fprintf('  Years where |delta fval| < 0.01:  %d / %d\n', sum(abs(fval_fc - fval_em) < 0.01), numYears);
fprintf('  Correlation of omegas (All-E):  %.4f\n', corr(omega_fc(4,:)', omega_em(4,:)'));
fprintf('  Correlation of omegas (High-E): %.4f\n', corr(omega_fc(3,:)', omega_em(3,:)'));
fprintf('  Correlation of omegas (High-U): %.4f\n', corr(omega_fc(2,:)', omega_em(2,:)'));
fprintf('  Correlation of omegas (High-N): %.4f\n', corr(omega_fc(1,:)', omega_em(1,:)'));
fprintf('  Correlation of omegas (All-N):  %.4f\n', corr(omega_fc(5,:)', omega_em(5,:)'));

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
