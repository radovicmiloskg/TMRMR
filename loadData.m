function [ DATA,LABELS ] = loadData( fileName )

load(fileName)

% For RSV group subject 020 had late symptoms and uninterpretable culture (Zaas et al.,2009). Excluding data:
if strcmp(fileName,'RSV.mat');
    newX(newX(:,2)==20,:)=[];
end

noOfGenes=newX(size(newX,1),1,1);
noOfTimePoints=size(newX,2)-3;
noOfPatient=size(newX,1)/noOfGenes;

DATA=zeros(noOfGenes,noOfTimePoints,noOfPatient);
LABELS=zeros(noOfPatient,1);

geneData=newX(:,4:size(newX,2));
geneData(geneData==0)=nan;

%% Interpolate missing values
for i=1:size(geneData,1)
    if any(isnan(geneData(i,:)))
        x=geneData(i,:);
        y=x;
        bd=isnan(x);
        gd=find(~bd);
        y(bd)=interp1(gd,x(gd),find(bd),'linear','extrap');
        geneData(i,:)=y;
    end
end

%% Set labels
for i=1:noOfPatient
    DATA(:,:,i)=geneData((i-1)*noOfGenes+1:i*noOfGenes,:);
    LABELS(i)=newX((i-1)*noOfGenes+1,3);
end
LABELS(LABELS<0)=0;

end

