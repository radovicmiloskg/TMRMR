function [geneID_TMRMR_C,geneID_TMRMR_M]=temporalMRMR(DATA,LABELS,genesToSelect)
%% RELEVANCE
alpha=0.3;
F=findCorrWithClassTemporal(DATA,LABELS);
[sortedF,indSortedF]=sort(F,'descend');

noOfGenes=size(DATA,1);
TopGenes=round(alpha*noOfGenes);
DATAreduced=DATA(indSortedF(1:TopGenes),:,:);
Fpool=[indSortedF(1:TopGenes) sortedF(1:TopGenes)];

%% REDUNCANCY
[C_cross_DTW,C_match_DTW]=findCorrAmongAttributesTemporal(DATAreduced);

%% TMRMR-C
Faccepted=zeros(genesToSelect,2);
Faccepted(1,:)=Fpool(1,:);
acceptedIndVec=zeros(1,TopGenes);
acceptedIndVec(1)=1;
acceptedIndVec=logical(acceptedIndVec);
TotalWc=0;
NoOfCorelations=0;
for i=2:genesToSelect
    Vf=1/i*(sum(Faccepted(:,2)) + Fpool(:,2));
    Vf(acceptedIndVec)=nan;
    
    NoOfCorelations=NoOfCorelations+(i-1);
    Wc=1/NoOfCorelations*(TotalWc + sum(abs(C_cross_DTW(acceptedIndVec,:)),1) )';
    Wc(acceptedIndVec)=nan;
    FCQ=Vf./Wc;
    
    [~,indFCQ]=nanmax(FCQ);    
    acceptedIndVec(indFCQ)=1;
    
    TotalWc=Wc(indFCQ)*NoOfCorelations;

    Faccepted(i,:)=Fpool(indFCQ,:);
end
geneID_TMRMR_C=Faccepted(:,1);
geneID_TMRMR_C(TopGenes+1:length(geneID_TMRMR_C))=nan;

%% TMRMR-M
Faccepted=zeros(genesToSelect,2);
Faccepted(1,:)=Fpool(1,:);
acceptedIndVec=zeros(1,TopGenes);
acceptedIndVec(1)=1;
acceptedIndVec=logical(acceptedIndVec);
TotalWc=0;
NoOfCorelations=0;
for i=2:genesToSelect
    Vf=1/i*(sum(Faccepted(:,2)) + Fpool(:,2));
    Vf(acceptedIndVec)=nan;    
 
    NoOfCorelations=NoOfCorelations+(i-1);
    Wc=1/NoOfCorelations*(TotalWc + sum(abs(C_match_DTW(acceptedIndVec,:)),1) )';
    Wc(acceptedIndVec)=nan;
    FCQ=Vf./Wc;
    
    [~,indFCQ]=nanmax(FCQ);   
    acceptedIndVec(indFCQ)=1;
    
    TotalWc=Wc(indFCQ)*NoOfCorelations;

    Faccepted(i,:)=Fpool(indFCQ,:);
end
geneID_TMRMR_M=Faccepted(:,1);
geneID_TMRMR_M(TopGenes+1:length(geneID_TMRMR_M))=nan;

end