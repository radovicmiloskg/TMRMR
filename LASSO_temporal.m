function [ LassoTemporal ] = LASSO_temporal( DATA,LABELS,genesToSelect )

addpath('.\MALSAR\MALSAR\functions\joint_feature_learning');
addpath('.\MALSAR\MALSAR\utils');
addpath('.\MALSAR\MALSAR\functions\Lasso');
addpath('.\MALSAR\MALSAR\functions\progression_model\TGL');

%% Rearange data:
noOfTimeSteps=size(DATA,2);

for i=1:noOfTimeSteps
    temp=DATA(:,i,:);
    temp=squeeze(temp);
    temp=temp';
    
    REARANGED_DATA{i}=temp;
    REARANGED_LABELS{i}=LABELS;
end

rho1=0;
rho2=0;
rho3=0.1;
currentSelection=[];

while sum(currentSelection)<genesToSelect
    [W, funcVal] = Least_TGL(REARANGED_DATA, REARANGED_LABELS, rho1, rho2, rho3);
    currentSelection=sum(W(:,1)~=0);
    rho3=0.5*rho3;
end

sumOfWeights=sum(abs(W),2);

[sortedW,ind]=sort(sumOfWeights,'descend');
LassoTemporal=ind(1:genesToSelect);

