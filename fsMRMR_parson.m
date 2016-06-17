function out = fsMRMR_parson( X, Y, param )
% function out = fsMRMR_parson(X, Y, param)
% 
% FCQ/FCD scheme according to MRMR
% X - the data, each row is an instance  
% Y - the class label in  form of 1 2 3 4 5 6 7 ...
% param.k - number of selected features
% param.pool - the number features will be considered in the second iteration 
% param.type = 1: FCQ, -1:FCD 

if nargin < 3
    param.k = 10; param.pool = 1000; param.type = 1;
end

if ~isfield(param, 'k')
    param.k = 10;
end

if ~isfield(param, 'pool')
    param.pool = 1000;
end

if ~isfield(param, 'type')
    param.type = 1;
end

%Number of features:
nF = size(X,2);

fea = zeros(param.k,1);

% correlation
t = abs(corr(X,Y));
t(isnan(t)) = 0;

% add features provide lowest redundancy
[orMI, idxs] = sort(t,'descend'); 
fea(1) = idxs(1); 
KMAX = min(param.pool,nF); 

idxleft = idxs(2:KMAX); 
mi_array = zeros(max(idxleft), param.k);

curlastfea = 1;

for ft = 2:param.k
   mi_array(idxleft,curlastfea) = abs(corr(X(:,idxleft),X(:,fea(curlastfea))));
   c_mi = mean(mi_array(idxleft, 1:ft-1),2);
   switch (param.type)
       case 1
           [tmp, tmpidx] = max(t(idxleft) ./ (c_mi + 0.00));
       case -1
           [tmp, tmpidx] = max(t(idxleft) - c_mi);
   end
   
   fea(ft) = idxleft(tmpidx); idxleft(tmpidx) = []; curlastfea = curlastfea + 1;
end

out.fList = fea; out.W = zeros(nF,1); out.W(fea) = param.k:-1:1; out.prf = -1;