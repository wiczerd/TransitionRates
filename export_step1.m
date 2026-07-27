% export_step1.m — Export Step 1 results to Excel
%
% Creates a clearly labeled Excel file with:
%   Sheet 1: Omegas (type shares) over time
%   Sheet 2: Ergodic distributions by type over time
%   Sheet 3: Unemployment, OLF, Employment rates by type
%   Sheet 4: Transition matrix parameters by type
%   Sheet 5: Convergence diagnostics

if ~exist('cfg','var')
    config_Paper2;
    cfg.mode = 'free_omega';
    cfg.Theta = 3; cfg.ThetaTot = 5;
    cfg.omegaEstimated = true; cfg.allE_estimated = false;
    cfg.outputPrefix = 'step1';
    cfg.outputDir = fullfile(cfg.dataDir, 'Step1_FreeOmega');
end

numYears = cfg.numYears;
years = cfg.yearsList(:);
Theta = cfg.Theta;
ThetaTot = cfg.ThetaTot;

%% Load all results
omega_ts = zeros(ThetaTot, numYears);
ergMat_ts = zeros(4, Theta, numYears);
pVecs_ts = zeros(12, Theta, numYears);
SSR2_ts = zeros(numYears, 1);
ef_ts = zeros(numYears, 1);

for yr = 1:numYears
    fname = fullfile(cfg.outputDir, sprintf('%s_%d.mat', cfg.outputPrefix, yr));
    if ~isfile(fname), warning('Missing: %s', fname); continue; end
    tmp = load(fname);
    r = tmp.res;
    omega_ts(:,yr) = r.omega;
    SSR2_ts(yr) = r.SSR2;
    ef_ts(yr) = r.exitflag;
    pVecs_ts(:,:,yr) = r.pVecs;
    for th = 1:Theta
        ergMat_ts(:,th,yr) = r.ergVec{th};
    end
end

xlsFile = fullfile(cfg.outputDir, 'Step1_FreeOmega_Results.xlsx');

%% Sheet 1: Omegas
typeNamesAll = {'All-N','High-N','High-U','High-E','All-E'};
% Reorder to match the 5-type omega vector:
% omega = [mover1=High-N, mover2=High-U, mover3=High-E, All-E, All-N]
% Display order: All-N, High-N, High-U, High-E, All-E
displayOrder = [5, 1, 2, 3, 4];

header_om = [{'Year'}, typeNamesAll];
data_om = [num2cell(years), num2cell(omega_ts(displayOrder,:)')];
writecell([header_om; data_om], xlsFile, 'Sheet', 'Omegas');

%% Sheet 2: Ergodic distributions
stateNames = {'N_OLF','U','E_short','E_long'};
typeNames = {'High-N','High-U','High-E'};
header_erg = {'Year'};
for th = 1:Theta
    for s = 1:4
        header_erg{end+1} = sprintf('%s_%s', typeNames{th}, stateNames{s});
    end
end
data_erg = zeros(numYears, Theta*4);
for yr = 1:numYears
    for th = 1:Theta
        data_erg(yr, (th-1)*4+1:th*4) = ergMat_ts(:,th,yr)';
    end
end
writecell([header_erg; num2cell([years, data_erg])], xlsFile, 'Sheet', 'Ergodic');

%% Sheet 3: Rates (U-rate, OLF-rate, E-rate by type + aggregate)
header_rates = {'Year'};
for th = 1:Theta
    header_rates{end+1} = sprintf('%s_Urate', typeNames{th});
end
header_rates{end+1} = 'Agg_Urate';
for th = 1:Theta
    header_rates{end+1} = sprintf('%s_OLFrate', typeNames{th});
end
header_rates{end+1} = 'Agg_OLFrate';
for th = 1:Theta
    header_rates{end+1} = sprintf('%s_Erate', typeNames{th});
end
header_rates{end+1} = 'Agg_Erate';

data_rates = zeros(numYears, numel(header_rates)-1);
for yr = 1:numYears
    om = omega_ts(:,yr);
    % Aggregate ergodic = sum over mover types weighted by omega + polar contributions
    ergAgg = zeros(4,1);
    for th = 1:Theta
        ergAgg = ergAgg + om(th) * ergMat_ts(:,th,yr);
    end
    ergAgg = ergAgg + [om(ThetaTot); 0; 0; om(ThetaTot-1)]; % All-N->state1, All-E->state4

    col = 0;
    % U-rate by type: erg(2) / (erg(2)+erg(3)+erg(4))
    for th = 1:Theta
        ev = ergMat_ts(:,th,yr);
        col = col+1;
        denom = ev(2)+ev(3)+ev(4);
        if denom > 0, data_rates(yr,col) = ev(2)/denom; end
    end
    % Aggregate U-rate
    col = col+1;
    denom = ergAgg(2)+ergAgg(3)+ergAgg(4);
    if denom > 0, data_rates(yr,col) = ergAgg(2)/denom; end

    % OLF-rate by type: erg(1)
    for th = 1:Theta
        col = col+1;
        data_rates(yr,col) = ergMat_ts(1,th,yr);
    end
    col = col+1;
    data_rates(yr,col) = ergAgg(1);

    % E-rate by type: erg(3)+erg(4)
    for th = 1:Theta
        col = col+1;
        data_rates(yr,col) = ergMat_ts(3,th,yr) + ergMat_ts(4,th,yr);
    end
    col = col+1;
    data_rates(yr,col) = ergAgg(3)+ergAgg(4);
end
writecell([header_rates; num2cell([years, data_rates])], xlsFile, 'Sheet', 'Rates');

%% Sheet 4: Transition parameters
% pVec = [p12,p13,p14, p22,p23,p24, p32,p33,p34, p42,p43,p44]
paramNames = {'p12','p13','p14','p22','p23','p24','p32','p33','p34','p42','p43','p44'};
header_params = {'Year'};
for th = 1:Theta
    for p = 1:12
        header_params{end+1} = sprintf('%s_%s', typeNames{th}, paramNames{p});
    end
end
data_params = zeros(numYears, Theta*12);
for yr = 1:numYears
    for th = 1:Theta
        data_params(yr, (th-1)*12+1:th*12) = pVecs_ts(:,th,yr)';
    end
end
writecell([header_params; num2cell([years, data_params])], xlsFile, 'Sheet', 'Parameters');

%% Sheet 5: Convergence
header_conv = {'Year','ExitFlag','Fval'};
data_conv = [years, ef_ts, SSR2_ts];
writecell([header_conv; num2cell(data_conv)], xlsFile, 'Sheet', 'Convergence');

fprintf('Excel saved: %s\n', xlsFile);
