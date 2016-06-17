function [optimal_K] = findOptimalK(INPUTS,LABELS,KFold)

K=[1 3 5 7];
indices = crossvalind('Kfold', LABELS, KFold);
KNN_ACC=zeros(length(K),1);

for i=1:length(K)
    for j=1:KFold
        testInd = (indices == j);
        trainInd = ~testInd;

        DATA_TRAIN=INPUTS;
        DATA_TRAIN(testInd,:)=[];

        LABELS_TRAIN=LABELS;
        LABELS_TRAIN(testInd)=[];

        DATA_TEST=INPUTS(testInd,:);
        LABELS_TEST=LABELS(testInd);

        KNN = fitcknn(DATA_TRAIN,LABELS_TRAIN,'NumNeighbors',K(i));
        predictedClass_KNN = predict(KNN,DATA_TEST);
        KNN_ACC(i)=KNN_ACC(i)+sum(predictedClass_KNN==LABELS_TEST);
    end
end

[~,ind]=max(KNN_ACC);

optimal_K=K(ind);

