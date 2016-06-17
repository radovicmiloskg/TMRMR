%% FRESH START
clear 
close all
clc
rng(0)

%% MEX DTW
mex dtw_c.c;

%% LOAD DATA
datasetName='H3N2'; % 'H3N2', 'Rhino' or 'RSV'
[DATA,LABELS]=loadData([datasetName '.mat']);

%% REDUCE NUMBER OF TIME POINTS
noOfTimeSteps=size(DATA,2);
timePointsToUse=noOfTimeSteps; % No reduction - noOfTimeSteps
if timePointsToUse<noOfTimeSteps
    step=noOfTimeSteps/(timePointsToUse-1);
    timePointsToKeep=[1];
    for i=1:timePointsToUse-2;
        timePointsToKeep=[timePointsToKeep round(timePointsToKeep(1)+i*step)];
    end
    timePointsToKeep=[timePointsToKeep noOfTimeSteps];
    DATA=DATA(:,timePointsToKeep,:);
end

%% SCALE DATA
% second arg: 'min-max' or 'z-score'
DATA=scaleData(DATA,'min-max');

%% SELECTED ATTRIBUTES STORAGE
FSS.MRMR=[];
FSS.F_statistic=[];
FSS.RELIEFF=[];
FSS.LASSO=[];
FSS.TMRMR_C=[];
FSS.TMRMR_M=[];

%% CROSS VALIDATION
NoOfClassifiers=3;
NoOfFeatureSelections=6;
genesToSelect=50;
genesForClassification=[1 10 20 30 40 50];
noOfPatient=size(DATA,3);

Kfold=5;
indices = crossvalind('Kfold', LABELS, Kfold);

RESULTS_ACC=zeros(length(genesForClassification),NoOfFeatureSelections,NoOfClassifiers);
RESULTS_SENS=zeros(length(genesForClassification),NoOfFeatureSelections,NoOfClassifiers);
RESULTS_SPEC=zeros(length(genesForClassification),NoOfFeatureSelections,NoOfClassifiers);

for i=1:Kfold
    testInd = (indices == i);
    trainInd = ~testInd;
    
    DATA_TRAIN=DATA;
    DATA_TRAIN(:,:,testInd)=[];
    
    LABELS_TRAIN=LABELS;
    LABELS_TRAIN(testInd)=[];
    
    DATA_TEST=DATA(:,:,testInd);
    LABELS_TEST=LABELS(testInd);
    
    % PREPROCESSING REQUIRED BY F-STATISTIC, MRMR and RELIEFF - DATA FLATTENING
    [DATA_TRAIN_FLATTEN,LABELS_TRAIN_FLATTEN]=flattenData(DATA_TRAIN,LABELS_TRAIN);
    
    % SELECT FEATURES - mRMR
    paramMRMR.k = genesToSelect; % Number of genes to select
    paramMRMR.type = 1; % FCQ    
    [geneID_MRMR] = fsMRMR_parson(DATA_TRAIN_FLATTEN,LABELS_TRAIN_FLATTEN, paramMRMR);
    geneID_MRMR=geneID_MRMR.fList;
    
    % SELECT FEATURES - F_statistic (ANOVA)
    geneID_F_statistic=calculateF_statistic(DATA_TRAIN_FLATTEN,LABELS_TRAIN_FLATTEN,genesToSelect);    
    
    % SELECT FEATURES - RELIEFF
    [geneID_RELIEFF,~] = relieff(DATA_TRAIN_FLATTEN,LABELS_TRAIN_FLATTEN,5,'method','classification');
    geneID_RELIEFF=geneID_RELIEFF(1:genesToSelect);

    % SELECT FEATURES - MLT-LASSO
    [geneID_LassoTemporal]=LASSO_temporal(DATA_TRAIN,LABELS_TRAIN,genesToSelect);
    
    % SELECT ATRIBUTES - TMRMR_C and TMRMR_M
    [geneID_TMRMR_C,geneID_TMRMR_M]=temporalMRMR(DATA_TRAIN,LABELS_TRAIN,genesToSelect);
    
    % CLASSIFICATION
    [currentACC,currentSENS,currentSPEC]=performClassification(DATA_TRAIN,LABELS_TRAIN,DATA_TEST,LABELS_TEST,genesForClassification,geneID_MRMR,geneID_F_statistic,geneID_RELIEFF,geneID_LassoTemporal,geneID_TMRMR_C,geneID_TMRMR_M);
    RESULTS_ACC=RESULTS_ACC+currentACC;
    RESULTS_SENS=RESULTS_SENS+currentSENS;
    RESULTS_SPEC=RESULTS_SPEC+currentSPEC;
    
    % SAVE SELECTED GENES FOR EACH FOLD
    FSS.MRMR=[FSS.MRMR geneID_MRMR];
    FSS.F_statistic=[FSS.F_statistic geneID_F_statistic];
    FSS.RELIEFF=[FSS.RELIEFF geneID_RELIEFF'];
    FSS.LASSO=[FSS.LASSO geneID_LassoTemporal];
    FSS.TMRMR_C=[FSS.TMRMR_C geneID_TMRMR_C];
    FSS.TMRMR_M=[FSS.TMRMR_M geneID_TMRMR_M];
end
RESULTS_ACC=RESULTS_ACC/noOfPatient;
RESULTS_SENS=RESULTS_SENS/sum(LABELS);
RESULTS_SPEC=RESULTS_SPEC/sum(~LABELS);

%% SAVE RESULTS
save(['RESULTS\' datasetName '_' num2str(timePointsToUse) '_ACC.mat'],'RESULTS_ACC')
save(['RESULTS\' datasetName '_' num2str(timePointsToUse) '_SENS.mat'],'RESULTS_SENS')
save(['RESULTS\' datasetName '_' num2str(timePointsToUse) '_SPEC.mat'],'RESULTS_SPEC')
save(['RESULTS\' datasetName '_' num2str(timePointsToUse) '_FSS.mat'],'FSS')

%% PLOT RESULTS
plotResults(RESULTS_ACC,genesForClassification,datasetName);
