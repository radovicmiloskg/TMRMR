function [] = plotResults( RESULTS,genesForClassification,str )

figure('units','normalized','outerposition',[0 0 1 1])
annotation('textbox',[0.49 0.9 0.1 0.1],'String',str,'FontSize',16,'LineStyle','none');

idx=isnan(RESULTS(:,1,1));
genesForClassification(idx)=[];
RESULTS(idx,:,:)=[];

subplot(3,1,1)
plot(genesForClassification,RESULTS(:,1,1),'b',genesForClassification,RESULTS(:,2,1),'y',genesForClassification,RESULTS(:,3,1),'r',genesForClassification,RESULTS(:,4,1),'m',genesForClassification,RESULTS(:,5,1),'g',genesForClassification,RESULTS(:,6,1),'k')
xlabel('No of genes')
ylabel('Accuracy');
legend('mRMR','F-statistic','RELIEFF','MT-LASSO','TMRMR-C','TMRMR-M');
title('KNN')
axis([0 genesForClassification(end) 0 1.1]);

subplot(3,1,2)
plot(genesForClassification,RESULTS(:,1,2),'b',genesForClassification,RESULTS(:,2,2),'y',genesForClassification,RESULTS(:,3,2),'r',genesForClassification,RESULTS(:,4,2),'m',genesForClassification,RESULTS(:,5,2),'g',genesForClassification,RESULTS(:,6,2),'k')
xlabel('No of genes')
ylabel('Accuracy');
legend('mRMR','F-statistic','RELIEFF','MT-LASSO','TMRMR-C','TMRMR-M');
title('NB')
axis([0 genesForClassification(end) 0 1.1]);

subplot(3,1,3)
plot(genesForClassification,RESULTS(:,1,3),'b',genesForClassification,RESULTS(:,2,3),'y',genesForClassification,RESULTS(:,3,3),'r',genesForClassification,RESULTS(:,4,3),'m',genesForClassification,RESULTS(:,5,3),'g',genesForClassification,RESULTS(:,6,3),'k')
xlabel('No of genes')
ylabel('Accuracy');
legend('mRMR','F-statistic','RELIEFF','MT-LASSO','TMRMR-C','TMRMR-M');
title('SVM')
axis([0 genesForClassification(end) 0 1.1]);

end

