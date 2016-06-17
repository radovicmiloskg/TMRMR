function [ C_cross_DTW,C_timeSteps_DTW ] = findCorrAmongAttributesTemporal( DATA )

noOfGenes=size(DATA,1);
noOfTimePoints=size(DATA,2);
noOfPatient=size(DATA,3);

C_cross_DTW=zeros(noOfGenes);
C_timeSteps_DTW=zeros(noOfGenes);

for i=1:noOfGenes
    Ai=DATA(i,:,:);
    Ai=reshape(Ai,[noOfTimePoints noOfPatient]);
    Ai=Ai';
    
    Ai_1=(Ai-repmat(nanmean(Ai,2),1,noOfTimePoints))./repmat(nanstd(Ai,0,2),1,noOfTimePoints);
    
    C_cross_DTW_Current=[];
    C_timeSteps_DTW_Current=[]; 

    parfor j=i+1:noOfGenes
        Aj=DATA(j,:,:);
        Aj=reshape(Aj,[noOfTimePoints noOfPatient]);
        Aj=Aj';
        
        Aj_1=(Aj-repmat(nanmean(Aj,2),1,noOfTimePoints))./repmat(nanstd(Aj,0,2),1,noOfTimePoints);
        
        w=50;
        DTW_cross=0;
        DTW_timeSteps=0;
        
        for k=1:noOfPatient
            DTW_timeSteps=DTW_timeSteps+dtw_c(Ai_1(k,~isnan(Ai_1(k,:)))',Aj_1(k,~isnan(Aj_1(k,:)))',w);
            for p=1:noOfPatient
                DTW_cross=DTW_cross+dtw_c(Ai_1(k,~isnan(Ai_1(k,:)))',Aj_1(p,~isnan(Aj_1(p,:)))',w);
            end
        end
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        DTW_cross=1/noOfPatient^2*DTW_cross;
        C_cross_DTW_Current(j)=1/DTW_cross;
        
        DTW_timeSteps=1/noOfPatient*DTW_timeSteps;
        C_timeSteps_DTW_Current(j)=1/DTW_timeSteps;     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    end
    
    C_cross_DTW(i+1:noOfGenes,i)=C_cross_DTW_Current(i+1:noOfGenes);
    C_cross_DTW(i,i+1:noOfGenes)=C_cross_DTW_Current(i+1:noOfGenes);
    
    C_timeSteps_DTW(i+1:noOfGenes,i)=C_timeSteps_DTW_Current(i+1:noOfGenes);
    C_timeSteps_DTW(i,i+1:noOfGenes)=C_timeSteps_DTW_Current(i+1:noOfGenes);

end