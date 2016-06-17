function [ F_statistic ] = calculateF_statistic( DATA,LABELS,genesToSelect )

noOfGenes=size(DATA,2);
F=zeros(noOfGenes,1);
classes=unique(LABELS);
m=length(unique(LABELS)); % Number of classes;

for i=1:noOfGenes
    Acurrent=DATA(:,i);
    AmeanAll=nanmean(Acurrent);
    AmeanClass=[];
    n_Ck=[];
    denom=0;
    
    for k=1:m
        AmeanClass(k)=nanmean(Acurrent(LABELS==classes(k)));
        n_Ck(k)=sum(LABELS==classes(k))-sum(isnan(Acurrent(LABELS==classes(k))));
        denom=denom+sum( (Acurrent(logical((LABELS==classes(k)) .* (~isnan(Acurrent))))-AmeanClass(k)).^2 );
    end
    
    temp(i)=abs(AmeanClass(1)-AmeanClass(2));
    
    N=sum(~isnan(Acurrent));
    denom=denom/(N-m);
    numer=n_Ck*((AmeanAll-AmeanClass).^2)'/(m-1);
    if denom==0
        F(i)=0;
    else
        F(i)=numer/denom;
    end
end

[~,indSortedF]=sort(F,'descend');
F_statistic=indSortedF(1:genesToSelect);

end

