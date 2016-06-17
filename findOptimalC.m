function [optimal_C] = findOptimalK(INPUTS,LABELS,KFold)

C=[10^(-3) 10^(-2) 10^(-1) 10^(0) 10^(1) 10^(2) 10^(3)];
indices = crossvalind('Kfold', LABELS, KFold);
SVM_ACC=zeros(length(C),1);

for i=1:length(C)
    for j=1:KFold
        testInd = (indices == j);
        trainInd = ~testInd;

        DATA_TRAIN=INPUTS;
        DATA_TRAIN(testInd,:)=[];

        LABELS_TRAIN=LABELS;
        LABELS_TRAIN(testInd)=[];

        DATA_TEST=INPUTS(testInd,:);
        LABELS_TEST=LABELS(testInd);
        
        SVM = fitcsvm(DATA_TRAIN,LABELS_TRAIN,'BoxConstraint',C(i));
        predictedClass_SVM = predict(SVM,DATA_TEST);
        SVM_ACC(i)=SVM_ACC(i)+sum(predictedClass_SVM==LABELS_TEST);
    end
end

[~,ind]=max(SVM_ACC);

optimal_C=C(ind);

