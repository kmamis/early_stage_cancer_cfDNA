clear all; close all;

%% -------- Load --------
load('cfDNA_Cohen_2018.mat');  % provides cfDNA_healthy and per-cancer vectors

%% -------- Vectors --------
H = cfDNA_healthy(:);

C = [cfDNA_lung(:);
    cfDNA_breast_I(:);
    cfDNA_breast_II_III(:);
    cfDNA_colorectal(:);
    cfDNA_pancreatic(:);
    cfDNA_ovarian(:);
    cfDNA_esophageal(:);
    cfDNA_stomach(:);
    cfDNA_liver(:)
];

%% -------- ksdensity (nonnegative support, reflection) --------
[f,  x ] = ksdensity(H, 'Support','nonnegative', 'BoundaryCorrection','reflection');
[ff, xx] = ksdensity(C, 'Support','nonnegative', 'BoundaryCorrection','reflection');

%% -------- Plot --------
figure; hold on
plot(x,  f,  'LineWidth', 2.5)
plot(xx, ff, 'LineWidth', 2.5)
xlim([0 50])

ax = gca;
ax.FontSize  = 20;
ax.LineWidth = 1.5;
ax.TickDir   = 'out';
box off

xlabel('cfDNA concentration (ng/mL)')
ylabel('probability distribution')

legend( ...
    sprintf('healthy individuals, {\\it N} = %d', numel(H)), ...
    sprintf('cancer patients, {\\it N} = %d', numel(C)), ...
    'FontSize', 16, 'Location','northeast' ...
);
