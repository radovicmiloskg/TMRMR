function [ DATA_FLATTEN,LABELS_FLATTEN ] = flattenData( DATA,LABELS )

noOfGenes=size(DATA,1);
noOfTimePoints=size(DATA,2);
noOfPatient=size(DATA,3);

DATA_FLATTEN=[];
LABELS_FLATTEN=[];
for i=1:noOfGenes
    Acurrent=DATA(i,:,:);
    Acurrent=reshape(Acurrent,[noOfTimePoints noOfPatient]);
    Acurrent=Acurrent';
    Acurrent=Acurrent(:);
    DATA_FLATTEN=[DATA_FLATTEN Acurrent];
end

LABELS_FLATTEN=repmat(LABELS,noOfTimePoints,1);

end

