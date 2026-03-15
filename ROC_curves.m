clear all

%% Parameters
k = 10;                               % Number of folds
nBoot = 1000;                         % For ROC confidence band and PPV/NPV CI
gridPoints = linspace(0,1,200);       % Specificity grid for interpolation

% Display names
cancerTypes  = {'Lung','Colorectal','Pancreatic','Ovarian','Esophageal','Stomach','Liver'};
breastStages = {'Breast_I','Breast_II_III'};
earlyStages  = {'Lung_I','Colorectal_I','Ovarian_I','Stomach_I'};

% Results table
results = table('Size',[numel(cancerTypes)+numel(breastStages)+numel(earlyStages) 17], ...
 'VariableTypes',{'string','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double'}, ...
 'VariableNames',{'CancerType','MeanAUC','AUC_CI_Lower','AUC_CI_Upper','OptCutoff','SensAtCutoff','Sens_CI_Lower','Sens_CI_Upper','SpecAtCutoff','Spec_CI_Lower','Spec_CI_Upper','MeanPPV','PPV_CI_Lower','PPV_CI_Upper','MeanNPV','NPV_CI_Lower','NPV_CI_Upper'});

%% Load data
S = load('cfDNA_Cohen_2018.mat');

% Choose healthy baseline 
healthy = S.cfDNA_healthy(:);

% Map display name -> variable in MAT
name2var = containers.Map( ...
    {'Lung','Colorectal','Pancreatic','Ovarian','Esophageal','Stomach','Liver'}, ...
    {'cfDNA_lung','cfDNA_colorectal','cfDNA_pancreatic','cfDNA_ovarian','cfDNA_esophageal','cfDNA_stomach','cfDNA_liver'} );

breast_map = containers.Map( ...
    {'Breast_I','Breast_II_III'}, ...
    {'cfDNA_breast_I','cfDNA_breast_II_III'} );

early_map = containers.Map( ...
    {'Lung_I','Colorectal_I','Ovarian_I','Stomach_I'}, ...
    {'cfDNA_lung_I','cfDNA_colorectal_I','cfDNA_ovarian_I','cfDNA_stomach_I'} );

%% Process all-stage cancers
row = 1;
for i = 1:numel(cancerTypes)
    varname = name2var(cancerTypes{i});
    cancer_vec = S.(varname)(:);
    results = processCancerGroup_vec(cancerTypes{i}, healthy, cancer_vec, k, nBoot, gridPoints, results, row);
    row = row + 1;
end

%% Process breast stages
for j = 1:numel(breastStages)
    varname = breast_map(breastStages{j});
    cancer_vec = S.(varname)(:);
    results = processCancerGroup_vec(breastStages{j}, healthy, cancer_vec, k, nBoot, gridPoints, results, row);
    row = row + 1;
end

%% Process Stage I cancers
for e = 1:numel(earlyStages)
    varname = early_map(earlyStages{e});
    cancer_vec = S.(varname)(:);
    results = processCancerGroup_vec(earlyStages{e}, healthy, cancer_vec, k, nBoot, gridPoints, results, row);
    row = row + 1;
end

%% Build the results table
fmtCI = @(m,l,u) sprintf('%.3f (%.3f, %.3f)', m, l, u);
n = height(results);

CancerType_col = strings(n,1);
AUC_col        = strings(n,1);
OptCutoff_col  = strings(n,1);
Sens_col       = strings(n,1);
Spec_col       = strings(n,1);
PPV_col        = strings(n,1);
NPV_col        = strings(n,1);

for i = 1:n
    CancerType_col(i) = string(results.CancerType(i));
    AUC_col(i)       = fmtCI(results.MeanAUC(i), results.AUC_CI_Lower(i), results.AUC_CI_Upper(i));
    OptCutoff_col(i) = sprintf('%.3f', results.OptCutoff(i));                % mean score threshold across folds
    Sens_col(i)      = fmtCI(results.SensAtCutoff(i), results.Sens_CI_Lower(i), results.Sens_CI_Upper(i));
    Spec_col(i)      = fmtCI(results.SpecAtCutoff(i), results.Spec_CI_Lower(i), results.Spec_CI_Upper(i));
    PPV_col(i)       = fmtCI(results.MeanPPV(i), results.PPV_CI_Lower(i), results.PPV_CI_Upper(i));
    NPV_col(i)       = fmtCI(results.MeanNPV(i), results.NPV_CI_Lower(i), results.NPV_CI_Upper(i));
end

DisplayResults = table( ...
    CancerType_col, AUC_col, OptCutoff_col, Sens_col, Spec_col, PPV_col, NPV_col, ...
    'VariableNames', { ...
        'Cancer Type', ...
        'AUC (95% CI)', ...
        'Optimal Cutoff', ...
        'Sensitivity at Cutoff (95% CI)', ...
        'Specificity at Cutoff (95% CI)', ...
        'PPV at Cutoff (95% CI)', ...
        'NPV at Cutoff (95% CI)'});

%% ---------------- Local functions ----------------
function results = processCancerGroup_vec(groupName, healthy_vec, cancer_vec, k, nBoot, gridPoints, results, rowIndex)
    % Build labels/scores
    healthy_vec = healthy_vec(:);
    cancer_vec  = cancer_vec(:);
    labels = [zeros(numel(healthy_vec),1); ones(numel(cancer_vec),1)];
    scores = [healthy_vec; cancer_vec];

    % Stratified CV when possible
    if min([sum(labels==0), sum(labels==1)]) < 10
        cv = cvpartition(numel(labels),'KFold',k);
    else
        cv = cvpartition(labels,'KFold',k);
    end

    aucVals = zeros(k,1); ppvVals = zeros(k,1); npvVals = zeros(k,1);
    sensVals = zeros(k,1); specVals = zeros(k,1);
    optThreshVals = zeros(k,1);                 % store per-fold score thresholds
    rocMatrix = zeros(k,numel(gridPoints));

    for fold = 1:k
        trainIdx = training(cv,fold); 
        testIdx  = test(cv,fold);

        [X,Y,T,AUC,optThresh,sens,spec,ppv,npv] = computeROC( ...
            scores(trainIdx), labels(trainIdx), scores(testIdx), labels(testIdx));

        aucVals(fold)=AUC; 
        ppvVals(fold)=ppv; 
        npvVals(fold)=npv;
        sensVals(fold)=sens; 
        specVals(fold)=spec;
        optThreshVals(fold) = optThresh;        % score cutoff chosen on training

        % Interpolate ROC of this fold on the common specificity grid
        specValsInterp = 1 - X; 
        [specValsInterp, uniqueIdx] = unique(specValsInterp,'stable');
        sensValsInterp = Y(uniqueIdx);
        rocMatrix(fold,:) = interp1(specValsInterp, sensValsInterp, gridPoints, 'linear','extrap');
    end

    % Confidence band across folds (mean curve + bootstrap band)
    [meanSens,lowerCI,upperCI] = computeConfidenceBand(rocMatrix,nBoot);

    % One displayed operating point from mean ROC (for Sens/Spec entries)
    Jmean = meanSens + gridPoints - 1; 
    [~, idxMaxJ] = max(Jmean);
    specAtCutoff = gridPoints(idxMaxJ);
    sensAtCutoff = meanSens(idxMaxJ);

    % Summary across folds
    meanAUC = mean(aucVals); 
    ciAUC   = prctile(aucVals,[2.5 97.5]);

    % PPV/NPV bootstrap across folds
    bootPPV = zeros(nBoot,1); 
    bootNPV = zeros(nBoot,1);
    for b = 1:nBoot
        idx = randsample(k,k,true);
        bootPPV(b) = mean(ppvVals(idx));
        bootNPV(b) = mean(npvVals(idx));
    end
    ciPPV = prctile(bootPPV,[2.5 97.5]); 
    ciNPV = prctile(bootNPV,[2.5 97.5]);

    % Sens/Spec CI across folds
    [ciSens, ciSpec] = computeSensSpecCI(sensVals, specVals, nBoot);

    % Optimal Cutoff reported = MEAN score threshold across folds
    optCutoff = mean(optThreshVals);

    % Fill results by row 
    results(rowIndex,:) = {groupName, ...
        meanAUC, ciAUC(1), ciAUC(2), ...
        optCutoff, ...
        sensAtCutoff, ciSens(1), ciSens(2), ...
        specAtCutoff, ciSpec(1), ciSpec(2), ...
        mean(ppvVals), ciPPV(1), ciPPV(2), ...
        mean(npvVals), ciNPV(1), ciNPV(2)};

    % --------- Display ROC figure ---------
    col = hex2rgb(color_for(groupName));
    fig = figure('Visible','on','Units','inches','Position',[1 1 6.93 5.2]);  % visible window
    xFill = [100*(1-gridPoints) 100*(1-fliplr(gridPoints))];
    yFill = [100*lowerCI       100*fliplr(upperCI)];
    fill(xFill,yFill,col,'FaceAlpha',0.15,'EdgeColor',col); hold on;
    plot(100*(1-gridPoints),100*meanSens,'Color',col,'LineWidth',2);
    plot([0 100],[0 100],'k','LineWidth',1.5);
    plot(100*(1-specAtCutoff),100*sensAtCutoff,'ro','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','r');
    set(gca,'FontSize',20,'LineWidth',2.5,'TickDir','out','Box','off');
    axis([0 100 0 100]);
    set(gca,'XTick',0:20:100,'XTickLabel',flip(0:20:100));
end

function [X,Y,T,AUC,optThresh,sens,spec,ppv,npv] = computeROC(trainScores,trainLabels,testScores,testLabels)
    % Build ROC on training and pick Youden-J threshold (score domain)
    [X,Y,T,AUC] = perfcurve(trainLabels,trainScores,1);
    sensitivity = Y; specificity = 1 - X;
    J = sensitivity + specificity - 1;
    [~,idx]  = max(J);
    optThresh = T(idx);                      % score threshold on training

    % Apply to test data
    predicted = testScores >= optThresh;
    TP = sum(predicted==1 & testLabels==1);
    TN = sum(predicted==0 & testLabels==0);
    FP = sum(predicted==1 & testLabels==0);
    FN = sum(predicted==0 & testLabels==1);
    sens = TP/(TP+FN);
    spec = TN/(TN+FP);
    ppv  = TP/(TP+FP);
    npv  = TN/(TN+FN);
end

function [meanSens,lowerCI,upperCI] = computeConfidenceBand(rocMatrix,nBoot)
    nFolds = size(rocMatrix,1);
    bootSens = zeros(nBoot,size(rocMatrix,2));
    for b = 1:nBoot
        idx = randsample(nFolds,nFolds,true);
        bootSens(b,:) = mean(rocMatrix(idx,:),1);
    end
    meanSens = mean(bootSens,1);
    lowerCI  = prctile(bootSens,2.5);
    upperCI  = prctile(bootSens,97.5);
end

function [ciSens, ciSpec] = computeSensSpecCI(sensVals, specVals, nBoot)
    n = numel(sensVals);
    bootSensCutoff = zeros(nBoot,1);
    bootSpecCutoff = zeros(nBoot,1);
    for b = 1:nBoot
        idx = randsample(n, n, true);
        bootSensCutoff(b) = mean(sensVals(idx));
        bootSpecCutoff(b) = mean(specVals(idx));
    end
    ciSens = prctile(bootSensCutoff,[2.5 97.5]);
    ciSpec = prctile(bootSpecCutoff,[2.5 97.5]);
end

% --------- Helper: map group name to color (hex string) ---------
function hex = color_for(name)
    m = containers.Map( ...
      {'Lung','Colorectal','Pancreatic','Ovarian','Esophageal','Stomach','Liver', ...
       'Breast_I','Breast_II_III', ...
       'Lung_I','Colorectal_I','Ovarian_I','Stomach_I'}, ...
      {'#FFC107','#bd0013','#bd00a2','#72bd00','#7E2F8E','#2f8e7e','#bd4b00', ...
       '#7e391b','#fd7336', ...
       '#FFC107','#bd0013','#72bd00','#2f8e7e'} );
    hex = m(name);
end

% --------- Helper: hex color to RGB triplet ---------
function rgb = hex2rgb(hex)
    hex = char(hex);
    if hex(1) == '#', hex(1) = []; end
    rgb = [hex2dec(hex(1:2)) hex2dec(hex(3:4)) hex2dec(hex(5:6))] / 255;
end
