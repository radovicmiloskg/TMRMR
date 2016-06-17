function [F]=findCorrWithClassTemporal(DATA,LABELS)

noOfGenes=size(DATA,1);
noOfTimePoints=size(DATA,2);
noOfPatient=size(DATA,3);

F=zeros(noOfGenes,noOfTimePoints);

classes=unique(LABELS);
m=length(unique(LABELS)); % Number of classes;

for i=1:noOfGenes
    Acurrent=DATA(i,:,:);
    Acurrent=reshape(Acurrent,[noOfTimePoints noOfPatient]);
    Acurrent=Acurrent';
    
    for j=1:noOfTimePoints
        AmeanAll=nanmean(Acurrent(:,j));
        AmeanClass=[];
        n_Ck=[];
        denom=0;
        for k=1:m
            AmeanClass(k)=nanmean(Acurrent(LABELS==classes(k),j));
            n_Ck(k)=sum(LABELS==classes(k))-sum(isnan(Acurrent(LABELS==classes(k),j)));
            
            denom=denom+sum((Acurrent( logical((LABELS==classes(k)) .* (~isnan(Acurrent(:,j)))) ,j)-AmeanClass(k)).^2);
        end
        N=sum(~isnan(Acurrent(:,j)));
        denom=denom/(N-m);
        numer=n_Ck*((AmeanAll-AmeanClass).^2)'/(m-1);
        if denom==0
            F(i,j)=0;
        else
            F(i,j)=numer/denom;
        end
    end
end
F=1/noOfTimePoints*sum(F,2);
end

