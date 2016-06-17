function [ DATA ] = scaleData( DATA,type )

noOfGenes=size(DATA,1);
noOfTimePoints=size(DATA,2);
noOfPatient=size(DATA,3);

%% MIN-MAX SCALING
if strcmp(type,'min-max')
    for i=1:noOfGenes
        currentGene=DATA(i,:,:);
        currentGene=currentGene(:);
        minCurrentGene=min(currentGene);
        maxCurrentGene=max(currentGene);
        DATA(i,:,:)=(DATA(i,:,:)-minCurrentGene)/(maxCurrentGene-minCurrentGene);
    end
end

%% Z-score SCALING
if strcmp(type,'z-score')
    for i=1:noOfGenes
        currentGene=DATA(i,:,:);
        currentGene=currentGene(:);
        currentGene=zscore(currentGene);
        currentGene=reshape(currentGene,[1 noOfTimePoints noOfPatient]);
        DATA(i,:,:)=currentGene;
    end
end

end

