clear all

%% ---------- Healthy calibration (from Materials & Methods) ----------
k    = 34.2;    % saturation scale (ng/mL)
v    = 1;
c_nl = 0.4;     % support cap for normalized S_nl
c_l  = 0.6;     % support cap for normalized S_l

% Healthy normalized shedding means (sum ~0.10)
meanS_nl_H = 0.02;
meanS_l_H  = 0.08;

% Healthy shedding CVs
CVS_nl = 0.68;
CVS_l  = 0.61;

% Monte-Carlo draws
Nmc = 5e5;

%% ---------- Load data ----------
% each matrix: col1 = C_tot, col2 = C_l ; we derive C_nl = C_tot - C_l
load('cfDNA_Mattox_2023.mat', 'healthy_no_outliers', 'ovarian_no_outliers', 'pancreatic');

% choose deconvolution method
deconv_method=2; %for Sun et al. QP; deconv_method=3 for Sun et al. NNLS; deconv_method=4 for Moss et al. QP; deconv_method=5 for Moss et al. NNLS

parseCohort = @(M) deal(M(:,1), M(:,deconv_method), M(:,1)-M(:,deconv_method)); % (C_tot, C_l, C_nl)

[Ctot_H, C_l_H, C_nl_H] = parseCohort(healthy_no_outliers);
[Ctot_OV, C_l_OV, C_nl_OV] = parseCohort(ovarian_no_outliers);
[Ctot_PA, C_l_PA, C_nl_PA] = parseCohort(pancreatic);

%% ---------- Define the two cancer types and run sequentially ----------
% pol1 = non-leukocyte multiplicative factor; pol2 = leukocyte multiplicative factor
configs = {
    struct('name','ovarian',    'pol1',3.0, 'pol2',2.2, 'Ctot',Ctot_OV, 'C_l',C_l_OV, 'C_nl',C_nl_OV), ...
    struct('name','pancreatic', 'pol1',5.0, 'pol2',2.5, 'Ctot',Ctot_PA, 'C_l',C_l_PA, 'C_nl',C_nl_PA)
};

% storage for summary
results = struct();

for it = 1:numel(configs)
    cfg = configs{it};
    cancerType = cfg.name;

    % Apply multiplicative increases to normalized shedding means
    meanS_nl = min(cfg.pol1 * meanS_nl_H, c_nl - 1e-9);
    meanS_l  = min(cfg.pol2 * meanS_l_H,  c_l  - 1e-9);

    % Rescaled-Beta parameters on [0,c] given mean m and CV
    betaParams = @(m,c,cv) deal( ...
        (m/c) * ( (m/c*(1 - m/c)) / ((cv*m)^2 / (c^2)) - 1 ), ...
        (1 - m/c) * ( (m/c*(1 - m/c)) / ((cv*m)^2 / (c^2)) - 1 ) );

    [alpha_nl, beta_nl] = betaParams(meanS_nl, c_nl, CVS_nl);
    [alpha_l,  beta_l ] = betaParams(meanS_l,  c_l,  CVS_l );

    if any([alpha_nl, beta_nl, alpha_l, beta_l] < 1)
        warning('[%s] Infeasible Beta params. Try smaller multipliers.', cancerType);
        continue
    end

    % Monte-Carlo sampling for normalized shedding
    S_nl = c_nl * betarnd(alpha_nl, beta_nl, 1, Nmc);
    S_l  = c_l  * betarnd(alpha_l,  beta_l,  1, Nmc);

    % Saturation transform to concentrations (feasible region S_nl+S_l<1)
    den  = v - (S_nl + S_l);
    mask = den > 0;
    C_nl = k * S_nl(mask) ./ den(mask);
    C_l  = k * S_l (mask) ./ den(mask);
    C_tot = C_nl + C_l;

    rho_model = corr(C_nl', C_l');

    % Analytical joint pdf for (C_nl, C_l)
    fun_joint = @(x,y) (k*v^2./(k + x + y).^3) .* ...
                       (1/c_nl).*betapdf(v*x./(c_nl*(k + x + y)), alpha_nl, beta_nl) .* ...
                       (1/c_l ).*betapdf(v*y./(c_l *(k + x + y)), alpha_l,  beta_l );

    % Evaluate PDFs on grid
    zz = (0:0.01:200)';
    f_l      = arrayfun(@(z) integral(@(x) fun_joint(x, z), 0, inf, 'RelTol',1e-6,'AbsTol',1e-9), zz);
    f_nonl   = arrayfun(@(z) integral(@(y) fun_joint(z, y), 0, inf, 'RelTol',1e-6,'AbsTol',1e-9), zz);
    f_tot    = arrayfun(@(z) integral(@(x) fun_joint(x, z - x), 0, z,   'RelTol',1e-6,'AbsTol',1e-9), zz);

    % CDFs for K-S tests
    FF_nonl = [zz, cumtrapz(zz, f_nonl)];
    FF_l    = [zz, cumtrapz(zz, f_l)];
    FF_tot  = [zz, cumtrapz(zz, f_tot)];

    % Cohort vectors (first two cols only used, non-l derived)
    Ctot_C = cfg.Ctot;  C_l_C = cfg.C_l;  C_nl_C = cfg.C_nl;

    % K-S tests (model vs cohort)
    [~, p_nonl] = kstest(C_nl_C, 'CDF', FF_nonl);
    [~, p_l    ] = kstest(C_l_C,  'CDF', FF_l);
    [~, p_tot  ] = kstest(Ctot_C, 'CDF', FF_tot);

    % Store
    results.(cancerType).pol1 = cfg.pol1;
    results.(cancerType).pol2 = cfg.pol2;
    results.(cancerType).p_nonl = p_nonl;
    results.(cancerType).p_l    = p_l;
    results.(cancerType).p_tot  = p_tot;
    results.(cancerType).rho_model = rho_model;
    results.(cancerType).zz = zz;
    results.(cancerType).f_nonl = f_nonl;
    results.(cancerType).f_l    = f_l;
    results.(cancerType).f_tot  = f_tot;

    % --------- Plots ---------
    % Non-leukocyte
    figure; hold on; box off; set(gca,'FontSize',22,'LineWidth',2,'TickDir','out');
    histogram(C_nl_C, 'Normalization','pdf','FaceColor',[0.85 0.35 0.10],'EdgeColor','none');
    plot(zz, f_nonl, 'Color',[0.49 0.18 0.56],'LineWidth',2.5);
    xlabel('non-leukocyte cfDNA (ng/mL)'); ylabel('pdf'); title(sprintf('%s: non-l', cancerType), 'Interpreter','none');

    % Leukocyte
    figure; hold on; box off; set(gca,'FontSize',22,'LineWidth',2,'TickDir','out');
    histogram(C_l_C, 'Normalization','pdf','FaceColor',[0.00 0.45 0.74],'EdgeColor','none');
    plot(zz, f_l, 'Color',[0.49 0.18 0.56],'LineWidth',2.5);
    xlabel('leukocyte cfDNA (ng/mL)'); ylabel('pdf'); title(sprintf('%s: l', cancerType), 'Interpreter','none');

    % Total
    figure; hold on; box off; set(gca,'FontSize',22,'LineWidth',2,'TickDir','out');
    histogram(Ctot_C, 'Normalization','pdf','FaceColor',[0.93 0.69 0.13],'EdgeColor','none');
    plot(zz, f_tot, 'Color',[0.49 0.18 0.56],'LineWidth',2.5);
    xlabel('total cfDNA (ng/mL)'); ylabel('pdf'); title(sprintf('%s: total', cancerType), 'Interpreter','none');

    % Console summary
    fprintf('[%s]  pol_nonl=%.2f  pol_l=%.2f  |  KS p_nonl=%.3g  p_l=%.3g  p_tot=%.3g  |  rho_model=%.2f\n', ...
            cancerType, cfg.pol1, cfg.pol2, p_nonl, p_l, p_tot, rho_model);
end

% results struct holds fields: results.ovarian and results.pancreatic