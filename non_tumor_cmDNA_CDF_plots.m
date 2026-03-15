clear all; close all;

%% ---------- Load data (cmDNA) ----------
%   cmDNA_healthy
%   cmDNA_lung_0, cmDNA_breast_0, cmDNA_colorectal_0, cmDNA_pancreatic_0,
%   cmDNA_esophageal_0, cmDNA_stomach_0, cmDNA_liver_0
load('cmDNA_Cohen_2018.mat');

%% ---------- Healthy baseline ------------
H = cmDNA_healthy(:);

%% ---------- α values and colors (non‑tumor cmDNA) ----------
alpha_all = struct( ...
    'lung',        1.92, ...
    'breast',      2.42, ...
    'colorectal',  2.36, ...
    'pancreatic',  4.01, ...
    'esophageal',  3.31, ...
    'stomach',     4.38, ...
    'liver',      12.37 ...
);

C = struct();
C.healthy    = '#0072BD';
C.lung       = '#EDB120';
C.breast     = '#fd7336';
C.colorectal = '#bd0013';
C.pancreatic = '#bd00a2';
C.esophageal = '#7E2F8E';
C.stomach    = '#2f8e7e';
C.liver      = '#bd4b00';

%% ---------- Panels ----------
panels = { ...
  'lung',       'cmDNA_lung_0'; ...
  'breast',     'cmDNA_breast_0'; ...
  'colorectal', 'cmDNA_colorectal_0'; ...
  'pancreatic', 'cmDNA_pancreatic_0'; ...
  'esophageal', 'cmDNA_esophageal_0'; ...
  'stomach',    'cmDNA_stomach_0'; ...
  'liver',      'cmDNA_liver_0' ...
};

%% ---------- Figure style ----------
XMIN = 1e-4; XMAX = 1e2;     % adjust if you prefer
YT   = 0:0.2:1;
XT   = [1e-4, 1e-2, 1e0, 1e2];
LW   = 5.5; FS = 35; AXLW = 4;

%% ---------- Draw ECDFs ----------
for i = 1:size(panels,1)
    label   = panels{i,1};
    varname = panels{i,2};
    if ~exist(varname,'var'), continue; end

    X0 = eval(varname);   % this is already the "non‑tumor" subset by your naming
    X0 = X0(:);

    figure
    ax = gca; hold on;
    ax.FontSize = FS; ax.LineWidth = AXLW; ax.TickDir = 'out'; box off;

    % healthy (blue)
    q = cdfplot(H);                       set(q,'LineWidth',LW,'Color',C.healthy);

    % α · healthy (solid, cancer color)
    a = alpha_all.(label);
    q = cdfplot(a * H);                   set(q,'LineWidth',LW,'Color',C.(label));

    % cancer non‑tumor cmDNA (dashed, same color)
    q = cdfplot(X0);                      set(q,'LineWidth',LW,'Color',C.(label),'LineStyle',':');

    set(gca,'XScale','log'); xlim([XMIN XMAX]); ylim([0 1]);
    yticks(YT); xticks(XT);
    xlabel(''); ylabel(''); title(strrep(label,'_','\_'));
    grid off;

end