function [ ACCURACY,SENSITIVITY,SPECIFICITY ] = performClassification( DATA_TRAIN,LABELS_TRAIN,DATA_TEST,LABELS_TEST,genesForClassification,geneID_MRMR,geneID_F_statistic,geneID_RELIEFF,geneID_LassoTemporal,geneID_TMRMR_C,geneID_TMRMR_M)
rng(0)
noOfPatient_train=size(DATA_TRAIN,3);
noOfPatient_test=size(DATA_TEST,3);
nestedKFold=5;

ind=0;
for k=genesForClassification
    %% TRAIN DATA PREPARING
    inputsMRMR_train=DATA_TRAIN(geneID_MRMR(1:k),:,:);
    inputsF_statistic_train=DATA_TRAIN(geneID_F_statistic(1:k),:,:);
    inputsRELIEFF_train=DATA_TRAIN(geneID_RELIEFF(1:k),:,:);
    inputsTemporal_Lasso_train=DATA_TRAIN(geneID_LassoTemporal(1:k),:,:); 
    inputsTemporal_cross_DWT_train=DATA_TRAIN(geneID_TMRMR_C(1:k),:,:);
    inputsTemporal_match_DWT_train=DATA_TRAIN(geneID_TMRMR_M(1:k),:,:);
    
    INPUTS_MRMR_train=[];
    INPUTS_F_statistic_train=[];
    INPUTS_RELIEFF_train=[];
    INPUTS_temporal_Lasso_train=[];      
    INPUTS_temporal_cross_DTW_train=[];
    INPUTS_temporal_match_DTW_train=[];  

    for i=1:noOfPatient_train
        tempMRMR=inputsMRMR_train(:,:,i);
        tempF_statistic=inputsF_statistic_train(:,:,i);
        tempRELIEFF=inputsRELIEFF_train(:,:,i);
        temp_Lasso=inputsTemporal_Lasso_train(:,:,i);     
        temp_cross_DWT=inputsTemporal_cross_DWT_train(:,:,i);
        temp_match_DWT=inputsTemporal_match_DWT_train(:,:,i); 
        
        INPUTS_MRMR_train=[INPUTS_MRMR_train;tempMRMR(:)'];
        INPUTS_F_statistic_train=[INPUTS_F_statistic_train;tempF_statistic(:)'];
        INPUTS_RELIEFF_train=[INPUTS_RELIEFF_train;tempRELIEFF(:)'];            
        INPUTS_temporal_Lasso_train=[INPUTS_temporal_Lasso_train;temp_Lasso(:)'];    
        INPUTS_temporal_cross_DTW_train=[INPUTS_temporal_cross_DTW_train;temp_cross_DWT(:)'];
        INPUTS_temporal_match_DTW_train=[INPUTS_temporal_match_DTW_train;temp_match_DWT(:)'];
    end
    
    %% TEST DATA PREPARING    
    inputsMRMR_test=DATA_TEST(geneID_MRMR(1:k),:,:);
    inputsF_statistic_test=DATA_TEST(geneID_F_statistic(1:k),:,:);
    inputsRELIEFF_test=DATA_TEST(geneID_RELIEFF(1:k),:,:);
    inputsTemporal_Lasso_test=DATA_TEST(geneID_LassoTemporal(1:k),:,:);     
    inputsTemporal_cross_DWT_test=DATA_TEST(geneID_TMRMR_C(1:k),:,:);
    inputsTemporal_match_DWT_test=DATA_TEST(geneID_TMRMR_M(1:k),:,:); 
  
    INPUTS_MRMR_test=[];
    INPUTS_F_statistic_test=[];
    INPUTS_RELIEFF_test=[];
    INPUTS_temporal_Lasso_test=[];
    INPUTS_temporal_cross_DTW_test=[];
    INPUTS_temporal_match_DTW_test=[];

    for i=1:noOfPatient_test
        tempMRMR=inputsMRMR_test(:,:,i);
        tempF_statistic=inputsF_statistic_test(:,:,i);
        tempRELIEFF=inputsRELIEFF_test(:,:,i);
        temp_Lasso=inputsTemporal_Lasso_test(:,:,i);  
        temp_cross_DWT=inputsTemporal_cross_DWT_test(:,:,i);
        temp_match_DWT=inputsTemporal_match_DWT_test(:,:,i);  

        INPUTS_MRMR_test=[INPUTS_MRMR_test;tempMRMR(:)'];
        INPUTS_F_statistic_test=[INPUTS_F_statistic_test;tempF_statistic(:)'];
        INPUTS_RELIEFF_test=[INPUTS_RELIEFF_test;tempRELIEFF(:)'];  
        INPUTS_temporal_Lasso_test=[INPUTS_temporal_Lasso_test;temp_Lasso(:)'];  
        INPUTS_temporal_cross_DTW_test=[INPUTS_temporal_cross_DTW_test;temp_cross_DWT(:)'];
        INPUTS_temporal_match_DTW_test=[INPUTS_temporal_match_DTW_test;temp_match_DWT(:)']; 
    end
   
    %% KNN    
    % mRMR
    [optimal_K]=findOptimalK(INPUTS_MRMR_train,LABELS_TRAIN,nestedKFold);
    KNN_MRMR = fitcknn(INPUTS_MRMR_train,LABELS_TRAIN,'NumNeighbors',optimal_K);
    predictedClass_KNN_MRMR = predict(KNN_MRMR,INPUTS_MRMR_test);
    KNN_MRMR_ACC=sum(predictedClass_KNN_MRMR==LABELS_TEST);
    KNN_MRMR_SENS=sum(predictedClass_KNN_MRMR(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    KNN_MRMR_SPEC=sum(predictedClass_KNN_MRMR(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));   
    % F-statistic
    [optimal_K]=findOptimalK(INPUTS_F_statistic_train,LABELS_TRAIN,nestedKFold);
    KNN_F_statistic = fitcknn(INPUTS_F_statistic_train,LABELS_TRAIN,'NumNeighbors',optimal_K);
    predictedClass_KNN_F_statistic = predict(KNN_F_statistic,INPUTS_F_statistic_test);
    KNN_F_statistic_ACC=sum(predictedClass_KNN_F_statistic==LABELS_TEST);
    KNN_F_statistic_SENS=sum(predictedClass_KNN_F_statistic(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    KNN_F_statistic_SPEC=sum(predictedClass_KNN_F_statistic(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));      
    % RELIEFF
    [optimal_K]=findOptimalK(INPUTS_RELIEFF_train,LABELS_TRAIN,nestedKFold);
    KNN_RELIEFF = fitcknn(INPUTS_RELIEFF_train,LABELS_TRAIN,'NumNeighbors',optimal_K);
    predictedClass_KNN_RELIEFF = predict(KNN_RELIEFF,INPUTS_RELIEFF_test);
    KNN_RELIEFF_ACC=sum(predictedClass_KNN_RELIEFF==LABELS_TEST);
    KNN_RELIEFF_SENS=sum(predictedClass_KNN_RELIEFF(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    KNN_RELIEFF_SPEC=sum(predictedClass_KNN_RELIEFF(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));   
    % MT-LASSO
    [optimal_K]=findOptimalK(INPUTS_temporal_Lasso_train,LABELS_TRAIN,nestedKFold);    
    KNN_Lasso = fitcknn(INPUTS_temporal_Lasso_train,LABELS_TRAIN,'NumNeighbors',optimal_K);
    predictedClass_KNN_Lasso = predict(KNN_Lasso,INPUTS_temporal_Lasso_test);
    KNN_Lasso_ACC=sum(predictedClass_KNN_Lasso==LABELS_TEST);
    KNN_Lasso_SENS=sum(predictedClass_KNN_Lasso(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    KNN_Lasso_SPEC=sum(predictedClass_KNN_Lasso(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));  
    % TMRMR-C
    [optimal_K]=findOptimalK(INPUTS_temporal_cross_DTW_train,LABELS_TRAIN,nestedKFold);       
    KNN_cross_DTW = fitcknn(INPUTS_temporal_cross_DTW_train,LABELS_TRAIN,'NumNeighbors',optimal_K);
    predictedClass_KNN_cross_DTW = predict(KNN_cross_DTW,INPUTS_temporal_cross_DTW_test);
    KNN_cross_DTW_ACC=sum(predictedClass_KNN_cross_DTW==LABELS_TEST);
    KNN_cross_DTW_SENS=sum(predictedClass_KNN_cross_DTW(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    KNN_cross_DTW_SPEC=sum(predictedClass_KNN_cross_DTW(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));
    % TMRMR-M
    [optimal_K]=findOptimalK(INPUTS_temporal_match_DTW_train,LABELS_TRAIN,nestedKFold);       
    KNN_match_DTW = fitcknn(INPUTS_temporal_match_DTW_train,LABELS_TRAIN,'NumNeighbors',optimal_K);
    predictedClass_KNN_match_DTW = predict(KNN_match_DTW,INPUTS_temporal_match_DTW_test);
    KNN_match_DTW_ACC=sum(predictedClass_KNN_match_DTW==LABELS_TEST);
    KNN_match_DTW_SENS=sum(predictedClass_KNN_match_DTW(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    KNN_match_DTW_SPEC=sum(predictedClass_KNN_match_DTW(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));  

    %% Naive-Bayes     
    % mRMR
    NB_MRMR = fitNaiveBayes(INPUTS_MRMR_train,LABELS_TRAIN);
    predictedClass_NB_MRMR = NB_MRMR.predict(INPUTS_MRMR_test);
    NB_MRMR_ACC=sum(predictedClass_NB_MRMR==LABELS_TEST);
    NB_MRMR_SENS=sum(predictedClass_NB_MRMR(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    NB_MRMR_SPEC=sum(predictedClass_NB_MRMR(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));   
    % F-statistic
    NB_F_statistic = fitNaiveBayes(INPUTS_F_statistic_train,LABELS_TRAIN);
    predictedClass_NB_F_statistic = NB_F_statistic.predict(INPUTS_F_statistic_test);
    NB_F_statistic_ACC=sum(predictedClass_NB_F_statistic==LABELS_TEST);
    NB_F_statistic_SENS=sum(predictedClass_NB_F_statistic(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    NB_F_statistic_SPEC=sum(predictedClass_NB_F_statistic(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));      
    % RELIEFF
    NB_RELIEFF = fitNaiveBayes(INPUTS_RELIEFF_train,LABELS_TRAIN);
    predictedClass_NB_RELIEFF = NB_RELIEFF.predict(INPUTS_RELIEFF_test);
    NB_RELIEFF_ACC=sum(predictedClass_NB_RELIEFF==LABELS_TEST);
    NB_RELIEFF_SENS=sum(predictedClass_NB_RELIEFF(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    NB_RELIEFF_SPEC=sum(predictedClass_NB_RELIEFF(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));   
    % MT-LASSO
    NB_Lasso = fitNaiveBayes(INPUTS_temporal_Lasso_train,LABELS_TRAIN);
    predictedClass_NB_Lasso = NB_Lasso.predict(INPUTS_temporal_Lasso_test);
    NB_Lasso_ACC=sum(predictedClass_NB_Lasso==LABELS_TEST);
    NB_Lasso_SENS=sum(predictedClass_NB_Lasso(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    NB_Lasso_SPEC=sum(predictedClass_NB_Lasso(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0)); 
    % TMRMR-C
    NB_cross_DTW = fitNaiveBayes(INPUTS_temporal_cross_DTW_train,LABELS_TRAIN);
    predictedClass_NB_cross_DTW = NB_cross_DTW.predict(INPUTS_temporal_cross_DTW_test);
    NB_cross_DTW_ACC=sum(predictedClass_NB_cross_DTW==LABELS_TEST);
    NB_cross_DTW_SENS=sum(predictedClass_NB_cross_DTW(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    NB_cross_DTW_SPEC=sum(predictedClass_NB_cross_DTW(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));
    % TMRMR-M
    NB_match_DTW = fitNaiveBayes(INPUTS_temporal_match_DTW_train,LABELS_TRAIN);
    predictedClass_NB_match_DTW = NB_match_DTW.predict(INPUTS_temporal_match_DTW_test);
    NB_match_DTW_ACC=sum(predictedClass_NB_match_DTW==LABELS_TEST);
    NB_match_DTW_SENS=sum(predictedClass_NB_match_DTW(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    NB_match_DTW_SPEC=sum(predictedClass_NB_match_DTW(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0)); 
    
    %% Support vector machine
    % mRMR
    [optimal_C]=findOptimalC(INPUTS_MRMR_train,LABELS_TRAIN,nestedKFold);      
    SVM_MRMR = fitcsvm(INPUTS_MRMR_train,LABELS_TRAIN,'BoxConstraint',optimal_C);
    predictedClass_SVM_MRMR = predict(SVM_MRMR,INPUTS_MRMR_test);
    SVM_MRMR_ACC=sum(predictedClass_SVM_MRMR==LABELS_TEST);
    SVM_MRMR_SENS=sum(predictedClass_SVM_MRMR(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    SVM_MRMR_SPEC=sum(predictedClass_SVM_MRMR(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));   
    % F-statistic
    [optimal_C]=findOptimalC(INPUTS_F_statistic_train,LABELS_TRAIN,nestedKFold); 
    SVM_F_statistic = fitcsvm(INPUTS_F_statistic_train,LABELS_TRAIN,'BoxConstraint',optimal_C);
    predictedClass_SVM_F_statistic = predict(SVM_F_statistic,INPUTS_F_statistic_test);
    SVM_F_statistic_ACC=sum(predictedClass_SVM_F_statistic==LABELS_TEST);
    SVM_F_statistic_SENS=sum(predictedClass_SVM_F_statistic(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    SVM_F_statistic_SPEC=sum(predictedClass_SVM_F_statistic(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));      
    % RELIEFF
    [optimal_C]=findOptimalC(INPUTS_RELIEFF_train,LABELS_TRAIN,nestedKFold); 
    SVM_RELIEFF = fitcsvm(INPUTS_RELIEFF_train,LABELS_TRAIN,'BoxConstraint',optimal_C);
    predictedClass_SVM_RELIEFF = predict(SVM_RELIEFF,INPUTS_RELIEFF_test);
    SVM_RELIEFF_ACC=sum(predictedClass_SVM_RELIEFF==LABELS_TEST);
    SVM_RELIEFF_SENS=sum(predictedClass_SVM_RELIEFF(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    SVM_RELIEFF_SPEC=sum(predictedClass_SVM_RELIEFF(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));
    % MT-LASSO
    [optimal_C]=findOptimalC(INPUTS_temporal_Lasso_train,LABELS_TRAIN,nestedKFold); 
    SVM_Lasso = fitcsvm(INPUTS_temporal_Lasso_train,LABELS_TRAIN,'BoxConstraint',optimal_C);
    predictedClass_SVM_Lasso = predict(SVM_Lasso,INPUTS_temporal_Lasso_test);
    SVM_Lasso_ACC=sum(predictedClass_SVM_Lasso==LABELS_TEST);
    SVM_Lasso_SENS=sum(predictedClass_SVM_Lasso(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    SVM_Lasso_SPEC=sum(predictedClass_SVM_Lasso(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0)); 
    % TMRMR-C
    [optimal_C]=findOptimalC(INPUTS_temporal_cross_DTW_train,LABELS_TRAIN,nestedKFold); 
    SVM_cross_DTW = fitcsvm(INPUTS_temporal_cross_DTW_train,LABELS_TRAIN,'BoxConstraint',optimal_C);
    predictedClass_SVM_cross_DTW = predict(SVM_cross_DTW,INPUTS_temporal_cross_DTW_test);
    SVM_cross_DTW_ACC=sum(predictedClass_SVM_cross_DTW==LABELS_TEST);
    SVM_cross_DTW_SENS=sum(predictedClass_SVM_cross_DTW(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    SVM_cross_DTW_SPEC=sum(predictedClass_SVM_cross_DTW(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0));
    % TMRMR-M
    [optimal_C]=findOptimalC(INPUTS_temporal_match_DTW_train,LABELS_TRAIN,nestedKFold); 
    SVM_match_DTW = fitcsvm(INPUTS_temporal_match_DTW_train,LABELS_TRAIN,'BoxConstraint',optimal_C);
    predictedClass_SVM_match_DTW = predict(SVM_match_DTW,INPUTS_temporal_match_DTW_test);
    SVM_match_DTW_ACC=sum(predictedClass_SVM_match_DTW==LABELS_TEST);
    SVM_match_DTW_SENS=sum(predictedClass_SVM_match_DTW(LABELS_TEST==1)==LABELS_TEST(LABELS_TEST==1));
    SVM_match_DTW_SPEC=sum(predictedClass_SVM_match_DTW(LABELS_TEST==0)==LABELS_TEST(LABELS_TEST==0)); 
    
    %% Calculate accuracy
    ind=ind+1;
    ACCURACY_KNN(ind,:)=[KNN_MRMR_ACC KNN_F_statistic_ACC KNN_RELIEFF_ACC KNN_Lasso_ACC KNN_cross_DTW_ACC KNN_match_DTW_ACC ];
    ACCURACY_NB(ind,:)=[NB_MRMR_ACC NB_F_statistic_ACC NB_RELIEFF_ACC NB_Lasso_ACC NB_cross_DTW_ACC NB_match_DTW_ACC];
    ACCURACY_SVM(ind,:)=[SVM_MRMR_ACC SVM_F_statistic_ACC SVM_RELIEFF_ACC SVM_Lasso_ACC SVM_cross_DTW_ACC SVM_match_DTW_ACC];
    
    SENSITIVITY_KNN(ind,:)=[KNN_MRMR_SENS KNN_F_statistic_SENS KNN_RELIEFF_SENS KNN_Lasso_SENS KNN_cross_DTW_SENS KNN_match_DTW_SENS ];
    SENSITIVITY_NB(ind,:)=[NB_MRMR_SENS NB_F_statistic_SENS NB_RELIEFF_SENS NB_Lasso_SENS NB_cross_DTW_SENS NB_match_DTW_SENS];
    SENSITIVITY_SVM(ind,:)=[SVM_MRMR_SENS SVM_F_statistic_SENS SVM_RELIEFF_SENS SVM_Lasso_SENS SVM_cross_DTW_SENS SVM_match_DTW_SENS];    
    
    SPECIFICITY_KNN(ind,:)=[KNN_MRMR_SPEC KNN_F_statistic_SPEC KNN_RELIEFF_SPEC KNN_Lasso_SPEC KNN_cross_DTW_SPEC KNN_match_DTW_SPEC ];
    SPECIFICITY_NB(ind,:)=[NB_MRMR_SPEC NB_F_statistic_SPEC NB_RELIEFF_SPEC NB_Lasso_SPEC NB_cross_DTW_SPEC NB_match_DTW_SPEC];
    SPECIFICITY_SVM(ind,:)=[SVM_MRMR_SPEC SVM_F_statistic_SPEC SVM_RELIEFF_SPEC SVM_Lasso_SPEC SVM_cross_DTW_SPEC SVM_match_DTW_SPEC]; 
    
end

ACCURACY(:,:,1)=ACCURACY_KNN;
ACCURACY(:,:,2)=ACCURACY_NB;
ACCURACY(:,:,3)=ACCURACY_SVM;

SENSITIVITY(:,:,1)=SENSITIVITY_KNN;
SENSITIVITY(:,:,2)=SENSITIVITY_NB;
SENSITIVITY(:,:,3)=SENSITIVITY_SVM;

SPECIFICITY(:,:,1)=SPECIFICITY_KNN;
SPECIFICITY(:,:,2)=SPECIFICITY_NB;
SPECIFICITY(:,:,3)=SPECIFICITY_SVM;

end