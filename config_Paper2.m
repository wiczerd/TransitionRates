% config_Paper2.m — Configuration for Paper 2 estimation pipeline
% All settings in one place. Run this script to populate cfg struct.

%% Data
cfg.dataDir      = fileparts(mfilename('fullpath'));  % same folder as this file
cfg.oldDataDir   = fullfile(fileparts(cfg.dataDir), ...
                   'Replication Files and ResultsTablesFigures Files 2 (1)');
cfg.inputFileNames = {'DataTimeSeriesYW.xlsx', 'DataTimeSeriesPW.xlsx', ...
                      'DataTimeSeriesYM.xlsx',  'DataTimeSeriesPM_1976_2023.xlsx'};
cfg.groupNum     = 4;            % 1=YW, 2=PW, 3=YM, 4=PM (primary: prime-age men)
cfg.dataRange    = 'J2:BD6562';  % columns J-BD = 47 years, rows 2-6562 = 6561 paths
cfg.yearStart    = 1976;
cfg.yearEnd      = 2024;
cfg.missingYears = [1977, 1993]; % no 8-month match possible
cfg.yearsList    = setdiff(cfg.yearStart:cfg.yearEnd, cfg.missingYears);
cfg.numYears     = numel(cfg.yearsList);  % 47
cfg.NumRespondents = [339039, 1255294, 344935, 1164770]; % from Weight_by_group.xlsx

%% Model specification
cfg.Nstates   = 4;              % hidden states per type
cfg.Nact      = 3;              % activities: E=1, U=2, N=3
cfg.Nmonth    = 8;              % months observed per sequence
cfg.Nstatepath = cfg.Nstates^cfg.Nmonth;  % 4^8 = 65536
cfg.Nactpath   = cfg.Nact^cfg.Nmonth;     % 3^8 = 6561
cfg.map        = [3, 2, 1, 1];  % state->activity: state1->N, state2->U, state3->E, state4->E
cfg.paramsPerType = 12;         % 3 free entries per row x 4 rows (p11 is residual)

%% Mode: 'free_omega' (Step 1) or 'fixed_omega' (Step 3)
cfg.mode = 'fixed_omega';       % <<< CHANGE THIS to switch steps

%% Sub-mode for fixed_omega: 'both_polar_hardwired' (3a) or 'allE_estimated' (3b)
cfg.fixedOmegaSubmode = 'allE_estimated';  % <<< CHANGE THIS for 3a vs 3b

% Derived type counts
switch cfg.mode
    case 'free_omega'
        cfg.Theta    = 3;       % mover types estimated
        cfg.ThetaTot = 5;       % total types (3 mover + All-N + All-E)
        cfg.omegaEstimated = true;
    case 'fixed_omega'
        switch cfg.fixedOmegaSubmode
            case 'both_polar_hardwired'
                cfg.Theta    = 3;   % 3 mover types, both polars hardwired
                cfg.ThetaTot = 5;
            case 'allE_estimated'
                cfg.Theta    = 4;   % 3 mover + All-E estimated
                cfg.ThetaTot = 5;   % + All-N hardwired
        end
        cfg.omegaEstimated = false;
end

%% All-E structural specification (only active when allE_estimated)
cfg.allE_estimated = strcmp(cfg.mode,'fixed_omega') && ...
                     strcmp(cfg.fixedOmegaSubmode,'allE_estimated');
% All-E is always the last estimated type (theta == cfg.Theta when allE_estimated)
% Structural zeros: state 3 (short-term E) removed from All-E
% In the 12-param vector [p12,p13,p14, p22,p23,p24, p32,p33,p34, p42,p43,p44]:
%   p13=0 (idx 2), p23=0 (idx 5), p32=p33=p34=0 (idx 7,8,9), p43=0 (idx 11)
cfg.allE_zeroIdx = [2, 5, 7, 8, 9, 11];  % indices within the 12-param block
cfg.allE_p44_lb  = 0.975;                 % lower bound on p44 (high persistence in long-term E)

%% Omega file for fixed-omega mode
cfg.projectRoot = fileparts(fileparts(cfg.dataDir));  % Topic_Heterogeneity_over_the_business_cycle
cfg.omegaFile = fullfile(cfg.projectRoot, ...
    '4 Model_FixedOmegas', 'Results 2025-09-23 4 types matrices', ...
    'V2 Nudge smoothed', 'Output files', 'omega_fixed.mat');

%% Warm-start file (may not exist — build_initial_guesses handles gracefully)
cfg.guessFile = fullfile(cfg.oldDataDir, 'guess_14.mat');

%% Old output directory (for comparison)
cfg.oldOutputDir = fullfile(cfg.projectRoot, ...
    '4 Model_FixedOmegas', 'Results 2025-09-23 4 types matrices', ...
    'V2 Nudge smoothed', 'Output files');

%% Optimization
cfg.MaxEvals = 20000;
cfg.MaxIter  = 10000;
cfg.OptTols  = 0.001;           % OptimalityTolerance
cfg.seed     = 100;             % random seed for initial guesses

%% Regularization (soft penalties encouraging type separation)
cfg.useRegularization = true;   % set false to disable
cfg.reg.lambda = 6e-3;          % per-observation weight
cfg.reg.margin = 5e-3;          % margin for dominance constraints
cfg.reg.tau    = 1e-3;          % softplus smoothing parameter

%% Bootstrap
cfg.Nrep = 1;                   % 1 = no bootstrap; >1 = bootstrap reps

%% Output
cfg.outputDir = cfg.dataDir;    % save per-year .mat files here
switch cfg.mode
    case 'free_omega'
        cfg.outputPrefix = 'output';
    case 'fixed_omega'
        switch cfg.fixedOmegaSubmode
            case 'both_polar_hardwired'
                cfg.outputPrefix = 'outputFO_3a';
            case 'allE_estimated'
                cfg.outputPrefix = 'outputFO_AllE';
        end
end
