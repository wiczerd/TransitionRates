% export_step1_orig.m — Export Step 1 results to Excel (from original code output)

if ~exist('outDir','var')
    rewriteDir = fileparts(mfilename('fullpath'));
    outDir = fullfile(rewriteDir, 'Step1_FreeOmega');
end
yearsList = setdiff(1976:2024, [1977 1993]);
numYears = numel(yearsList);
Theta = 3; ThetaTot = 5;

%% Load all results
omega_ts = zeros(ThetaTot, numYears);
ergMat_ts = zeros(4, Theta, numYears);
pVecs_ts = zeros(12, Theta, numYears);
SSR2_ts = zeros(numYears, 1);
ef_ts = zeros(numYears, 1);

for yr = 1:numYears
    fname = fullfile(outDir, sprintf('step1_%d.mat', yr));
    if ~isfile(fname), warning('Missing: %s', fname); continue; end
    tmp = load(fname);
    r = tmp.res;
    omega_ts(:,yr) = r.omega;
    SSR2_ts(yr) = r.SSR2;
    ef_ts(yr) = r.exitflag;
    pVecs_ts(:,:,yr) = r.pVecs;
    ergMat_ts(:,:,yr) = r.ergVecs;
end

xlsFile = fullfile(outDir, 'Step1_FreeOmega_Results.xlsx');

%% Sheet 1: Omegas
% omega layout after relabeling: [High-N, High-U, High-E, All-E, All-N]
% Display order: All-N, High-N, High-U, High-E, All-E
typeNamesAll = {'All-N','High-N','High-U','High-E','All-E'};
displayOrder = [5, 1, 2, 3, 4];
header_om = [{'Year'}, typeNamesAll];
data_om = [num2cell(yearsList(:)), num2cell(omega_ts(displayOrder,:)')];
writecell([header_om; data_om], xlsFile, 'Sheet', 'Omegas');

%% Sheet 2: Ergodic distributions
stateNames = {'N_OLF','U','E_short','E_long'};
typeNames = {'High-N','High-U','High-E'};
header_erg = {'Year'};
for th = 1:Theta
    for s = 1:4
        header_erg{end+1} = sprintf('%s_%s', typeNames{th}, stateNames{s}); %#ok<AGROW>
    end
end
data_erg = zeros(numYears, Theta*4);
for yr = 1:numYears
    for th = 1:Theta
        data_erg(yr, (th-1)*4+1:th*4) = ergMat_ts(:,th,yr)';
    end
end
writecell([header_erg; num2cell([yearsList(:), data_erg])], xlsFile, 'Sheet', 'Ergodic');

%% Sheet 3: Rates
header_rates = {'Year'};
for th = 1:Theta, header_rates{end+1} = sprintf('%s_Urate', typeNames{th}); end %#ok<AGROW>
header_rates{end+1} = 'Agg_Urate';
for th = 1:Theta, header_rates{end+1} = sprintf('%s_OLFrate', typeNames{th}); end %#ok<AGROW>
header_rates{end+1} = 'Agg_OLFrate';

data_rates = zeros(numYears, numel(header_rates)-1);
for yr = 1:numYears
    om = omega_ts(:,yr);
    ergAgg = zeros(4,1);
    for th = 1:Theta
        ergAgg = ergAgg + om(th) * ergMat_ts(:,th,yr);
    end
    ergAgg = ergAgg + [om(ThetaTot);0;0;om(ThetaTot-1)];
    col = 0;
    for th = 1:Theta
        ev = ergMat_ts(:,th,yr);
        col = col+1;
        denom = ev(2)+ev(3)+ev(4);
        if denom > 0, data_rates(yr,col) = ev(2)/denom; end
    end
    col = col+1;
    denom = ergAgg(2)+ergAgg(3)+ergAgg(4);
    if denom > 0, data_rates(yr,col) = ergAgg(2)/denom; end
    for th = 1:Theta
        col = col+1; data_rates(yr,col) = ergMat_ts(1,th,yr);
    end
    col = col+1; data_rates(yr,col) = ergAgg(1);
end
writecell([header_rates; num2cell([yearsList(:), data_rates])], xlsFile, 'Sheet', 'Rates');

%% Sheet 4: Convergence
header_conv = {'Year','ExitFlag','Fval'};
data_conv = [yearsList(:), ef_ts, SSR2_ts];
writecell([header_conv; num2cell(data_conv)], xlsFile, 'Sheet', 'Convergence');

fprintf('Excel saved: %s\n', xlsFile);
