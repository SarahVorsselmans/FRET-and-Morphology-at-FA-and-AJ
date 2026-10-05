function [pVal, starStr, testUsed] = runGroupStat(groupA, groupB, testType)
%RUNGROUPSTAT  Compare two groups and return p-value and significance stars.
%
%   [pVal, starStr, testUsed] = runGroupStat(groupA, groupB, testType)
%
% Inputs
%   groupA   : numeric vector (NaNs are removed internally)
%   groupB   : numeric vector (NaNs are removed internally)
%   testType : 'welch'       — Welch two-sample t-test  (default)
%              'mannwhitney' — Mann-Whitney U / Wilcoxon rank-sum test
%
% Outputs
%   pVal     : p-value (scalar); NaN if either group is empty
%   starStr  : significance string: '***' / '**' / '*' / 'ns' / ''
%   testUsed : char string describing the test actually run
%
% Significance thresholds
%   p < 0.001  →  ***
%   p < 0.01   →  **
%   p < 0.05   →  *
%   otherwise  →  ns (not significant — returned as empty string for clean plots)

if nargin < 3 || isempty(testType), testType = 'welch'; end

% Remove NaNs
a = groupA(~isnan(groupA));
b = groupB(~isnan(groupB));

% Guard: need at least 2 observations per group for a meaningful test
if numel(a) < 2 || numel(b) < 2
    pVal     = NaN;
    starStr  = '';
    testUsed = 'none (insufficient data)';
    return;
end

switch lower(string(testType))
    case 'mannwhitney'
        [pVal, ~] = ranksum(a, b);
        testUsed  = 'Mann-Whitney U (Wilcoxon rank-sum)';

    otherwise  % 'welch'
        [~, pVal] = ttest2(a, b, 'Vartype', 'unequal');
        testUsed  = 'Welch two-sample t-test';
end

% Map to significance stars
if     pVal < 0.001,  starStr = '***';
elseif pVal < 0.01,   starStr = '**';
elseif pVal < 0.05,   starStr = '*';
else,                 starStr = '';    % not significant — omit from plot
end
end
