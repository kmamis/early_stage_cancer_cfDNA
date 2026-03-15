clear all; close all;

%% ---------- Load data ----------
load('cfDNA_Cohen_2018.mat');  % contains: cfDNA_healthy, cfDNA for all cancer types, and stage-I subsets

%% ---------- Healthy baseline ------------
H = cfDNA_healthy(:);

%% ---------- α values and colors ----------
% All stages
alpha_all = struct( ...
    'lung',          1.29, ...
    'breast_I',      1.15, ...
    'breast_II_III', 1.79, ...
    'colorectal',    2.54, ...
    'pancreatic',    2.69, ...
    'ovarian',       4.30, ...
    'esophageal',    5.24, ...
    'stomach',       5.66, ...
    'liver',        12.63 ...
);

% Stage I
alpha_stageI = struct( ...
    'lung_I',       1.33, ...
    'colorectal_I', 2.39, ...
    'ovarian_I',    6.99, ...
    'stomach_I',    4.83 ...
);

% Colors
C = struct();
C.healthy        = '#0072BD';
C.lung           = '#EDB120';
C.breast_I       = '#fd7336';
C.breast_II_III  = '#7e391b';
C.colorectal     = '#bd0013';
C.pancreatic     = '#bd00a2';
C.ovarian        = '#72bd00';
C.esophageal     = '#7E2F8E';
C.stomach        = '#2f8e7e';
C.liver          = '#bd4b00';
C.lung_I         = C.lung;
C.colorectal_I   = C.colorectal;
C.ovarian_I      = C.ovarian;
C.stomach_I      = C.stomach;

%% ---------- Panels to draw ----------
panels_all = { ...
  'lung',          'cfDNA_lung'; ...
  'breast_I',      'cfDNA_breast_I'; ...
  'breast_II_III', 'cfDNA_breast_II_III'; ...
  'colorectal',    'cfDNA_colorectal'; ...
  'pancreatic',    'cfDNA_pancreatic'; ...
  'ovarian',       'cfDNA_ovarian'; ...
  'esophageal',    'cfDNA_esophageal'; ...
  'stomach',       'cfDNA_stomach'; ...
  'liver',         'cfDNA_liver' ...
};

panels_I = { ...
  'lung_I',       'cfDNA_lung_I'; ...
  'colorectal_I', 'cfDNA_colorectal_I'; ...
  'ovarian_I',    'cfDNA_ovarian_I'; ...
  'stomach_I',    'cfDNA_stomach_I' ...
};

%% ---------- Figure style ----------
XMIN = 1e-2; XMAX = 1e3;
YT   = 0:0.2:1;
XT   = [1e-2, 1e0, 1e2];
LW   = 5.5; FS = 35; AXLW = 4;

%% ---------- Helper to draw one panel ----------
draw_panel = @(lbl, varname, alpha_val) ...
    local_draw_ecdf_panel(lbl, varname, alpha_val, H, C, XMIN, XMAX, YT, XT, LW, FS, AXLW);

%% ---------- Draw ALL‑STAGES panels ----------
for i = 1:size(panels_all,1)
    label   = panels_all{i,1};
    varname = panels_all{i,2};
    if ~exist(varname,'var'), continue; end
    draw_panel(label, varname, alpha_all.(label));
end

%% ---------- Draw STAGE‑I panels ----------
for i = 1:size(panels_I,1)
    label   = panels_I{i,1};
    varname = panels_I{i,2};
    if ~exist(varname,'var'), continue; end
    draw_panel(label, varname, alpha_stageI.(label));
end

%% =======================================================================
%%                           Local function
%% =======================================================================
function local_draw_ecdf_panel(label, varname, alpha_val, H, C, XMIN, XMAX, YT, XT, LW, FS, AXLW)
    X = evalin('base', varname);  % pull vector from base workspace

    figure; ax = gca; hold on;
    ax.FontSize = FS; ax.LineWidth = AXLW; ax.TickDir = 'out'; box off;

    % Healthy (blue)
    q = cdfplot(H);            set(q,'LineWidth',LW,'Color',C.healthy);
    % α * Healthy (solid, color)
    clr = C.(label);           q = cdfplot(alpha_val * H); set(q,'LineWidth',LW,'Color',clr);
    % Cancer (dashed, same color)
    q = cdfplot(X);            set(q,'LineWidth',LW,'Color',clr,'LineStyle',':');

    set(gca,'XScale','log'); xlim([XMIN XMAX]); ylim([0 1]);
    yticks(YT); xticks(XT); xlabel(''); ylabel(''); title(strrep(label,'_','\_')); grid off;
end