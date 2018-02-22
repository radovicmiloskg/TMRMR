%% FRESH START
clear 
close all
clc
rng(0)

%% MEX DTW
mex dtw_c.c;

%% LOAD DATA
datasetName='GeneStudy001';
load('Dataset.mat'); %
% Dataset.mat should contain the following vatiables:
% 1) DATA - GxTxN matrix
% 2) LABELS - Nx1 vector
% where N is the number of subjects, G is the number of genes and T is the
% number of timesteps.

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

%% FEATURE SELECTION
genesToSelect=50;

% PREPROCESSING REQUIRED BY F-STATISTIC, MRMR and RELIEFF - DATA FLATTENING
[DATA_FLATTEN,LABELS_FLATTEN]=flattenData(DATA,LABELS);

% SELECT FEATURES - mRMR
paramMRMR.k = genesToSelect; % Number of genes to select
paramMRMR.type = 1; % FCQ    
[geneID_MRMR] = fsMRMR_parson(DATA_FLATTEN,LABELS_FLATTEN, paramMRMR);
FSS.MRMR=geneID_MRMR.fList;

% SELECT FEATURES - F_statistic (ANOVA)
FSS.F_statistic=calculateF_statistic(DATA_FLATTEN,LABELS_FLATTEN,genesToSelect);    

% SELECT FEATURES - RELIEFF
[geneID_RELIEFF,~] = relieff(DATA_FLATTEN,LABELS_FLATTEN,5,'method','classification');
FSS.RELIEFF=geneID_RELIEFF(1:genesToSelect);

% SELECT FEATURES - MLT-LASSO
[FSS.LASSO]=LASSO_temporal(DATA,LABELS,genesToSelect);

% SELECT ATRIBUTES - TMRMR_C and TMRMR_M
[FSS.TMRMR_C,FSS.TMRMR_M]=temporalMRMR(DATA,LABELS,genesToSelect);

%% SAVE SELECTED FEATURES
save(['RESULTS\' 'FSS_' datasetName '.mat'],'FSS');
