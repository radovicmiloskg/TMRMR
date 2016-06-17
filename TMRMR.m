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
FSS.TMRMR_C=[];
FSS.TMRMR_M=[];

%% FEATURE SELECTION
genesToSelect=50;
[FSS.TMRMR_C,FSS.TMRMR_M]=temporalMRMR(DATA,LABELS,genesToSelect);

%% SAVE SELECTED FEATURES
save(['RESULTS\' 'FSS_' datasetName '.mat'],'FSS');
