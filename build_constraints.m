function [A, B, Aeq, beq, lb, ub] = build_constraints(cfg)
% build_constraints  Build linear constraint matrices for fmincon.
%
%   [A, B, Aeq, beq, lb, ub] = build_constraints(cfg)
%
%   Parameter vector layout:
%     - Theta types x 12 transition params each = Theta*12 entries
%     - If free_omega mode: ThetaTot omega entries appended
%
%   Constraints:
%     A*x <= B:  each row of each transition matrix sums to <= 1
%                (with p_i1 = 1 - p_i2 - p_i3 - p_i4 >= 0)
%     Aeq*x = beq: omegas sum to 1 (free_omega mode only)
%     lb, ub: box bounds [0,1], with structural zeros for All-E

Theta    = cfg.Theta;
ThetaTot = cfg.ThetaTot;
ppt      = cfg.paramsPerType;  % 12
NParEst  = Theta * ppt;

if cfg.omegaEstimated
    Ntotal = NParEst + ThetaTot;
else
    Ntotal = NParEst;
end

% --- Inequality: A*x <= B (row sums of transition matrices) ---
Nrows = Theta * 4;  % 4 rows per transition matrix
A = zeros(Nrows, Ntotal);
B = ones(Nrows, 1);
for idx = 1:(Theta * 4)
    col_start = 3*(idx-1) + 1;
    A(idx, col_start:col_start+2) = 1;
end

% All-E: row 3 of its transition matrix is all zeros
if cfg.allE_estimated
    allE_row3 = (cfg.Theta - 1) * 4 + 3;  % row index for state 3 of last type
    B(allE_row3) = 0;
end

% --- Equality: Aeq*x = beq (omega sum to 1) ---
if cfg.omegaEstimated
    Aeq = [zeros(1, NParEst), ones(1, ThetaTot)];
    beq = 1;
else
    Aeq = [];
    beq = [];
end

% --- Box bounds ---
lb = zeros(Ntotal, 1);
ub = ones(Ntotal, 1);

% All-E structural zeros: fix certain params to 0 via ub = 0
if cfg.allE_estimated
    allE_base = (cfg.Theta - 1) * ppt;  % offset for the All-E type block
    for k = 1:numel(cfg.allE_zeroIdx)
        ub(allE_base + cfg.allE_zeroIdx(k)) = 0;
    end
    % p44 lower bound
    lb(allE_base + 12) = cfg.allE_p44_lb;
end

end
