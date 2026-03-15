clear all

%% Load dataset
load('cfDNA_Cohen_2018.mat');

% example
reference = cfDNA_colorectal(:);   % cancer group
target    = cfDNA_healthy(:);      % healthy group

%% Grid for multiplicative factor alpha
alpha_grid = 0:0.01:100;
N = numel(alpha_grid);

h = false(N,1);    % KS reject flag
p = zeros(N,1);    % p-values
D = zeros(N,1);    % KS distance between empirical cumulative distributions

%% Scan all alpha values 
for k = 1:N
    alpha = alpha_grid(k);

    % Compare alpha * healthy vs colorectal
    [h(k), p(k), D(k)] = kstest2(alpha * target, reference, 'Alpha', 0.05);
end

%% Best alpha = alpha that minimizes the KS distance
[~, idx_best] = min(D);
alpha_best = alpha_grid(idx_best);

%% KS 5% non-rejection interval
pass_idx = find(h == 0);            % h==0 -> KS DOES NOT reject
if isempty(pass_idx)
    alpha_lo = NaN;
    alpha_hi = NaN;
else
    alpha_lo = alpha_grid(pass_idx(1));
    alpha_hi = alpha_grid(pass_idx(end));
end

%% Display results (same as your original script)
disp('Best alpha (min KS distance):')
alpha_best

disp('Reject flag at alpha_best (0 = OK):')
h(idx_best)

disp('p-value at alpha_best:')
p(idx_best)

disp('5% KS non-rejection interval:')
[alpha_lo, alpha_hi]