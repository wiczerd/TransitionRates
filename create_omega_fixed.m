% create_omega_fixed.m — Create omega_fixed.mat for Phase 3 estimation
%
% Uses linear trend excluding 2020-2021 (option 3) from Phase 2 trend
% computation. Rearranges columns to match original code's expected format:
%   omega_fixed(year, :) = [calendarYear, All-N, High-N, High-U, High-E, All-E]

clc; clear;

rewriteDir = fileparts(mfilename('fullpath'));
load(fullfile(rewriteDir, 'Phase2_OmegaTrends', 'omega_trends.mat'), ...
    'yearsList', 'trend3_n', 'typeNames');

% trend3_n columns are [High-N, High-U, High-E, All-E, All-N] (from typeNames)
% Rearrange to original code format: [year, All-N, High-N, High-U, High-E, All-E]
omega_fixed = [yearsList(:), ...
    trend3_n(:,5), ...   % All-N
    trend3_n(:,1), ...   % High-N
    trend3_n(:,2), ...   % High-U
    trend3_n(:,3), ...   % High-E
    trend3_n(:,4)];      % All-E

save(fullfile(rewriteDir, 'omega_fixed.mat'), 'omega_fixed');
fprintf('Created omega_fixed.mat: %d years\n', size(omega_fixed,1));
fprintf('Columns: [year, All-N, High-N, High-U, High-E, All-E]\n');
fprintf('Trend method: Linear OLS excluding 2020-2021, normalized to sum=1\n\n');

% Print selected years
selYears = [1976 1980 1990 2000 2008 2010 2019 2020 2024];
fprintf('%6s %8s %8s %8s %8s %8s %8s\n', 'Year', 'All-N', 'High-N', 'High-U', 'High-E', 'All-E', 'Sum');
for i = 1:length(selYears)
    idx = find(omega_fixed(:,1) == selYears(i));
    if ~isempty(idx)
        row = omega_fixed(idx, 2:end);
        fprintf('%6d %8.4f %8.4f %8.4f %8.4f %8.4f %8.6f\n', ...
            omega_fixed(idx,1), row, sum(row));
    end
end
fprintf('\nSum check: min=%.8f, max=%.8f\n', ...
    min(sum(omega_fixed(:,2:end),2)), max(sum(omega_fixed(:,2:end),2)));
