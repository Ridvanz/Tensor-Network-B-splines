function yhat = evalspline(TN,tfeaturez,n,m)

% TN = trained tensor train

% Get basis vectors from test features
[ut] = basisvectors(tfeaturez,n,m);

% Contract all cores and basis vectors from left to right

[Nt, ~]=size(tfeaturez);
yhat = ones(Nt,1);

for i=1:size(TN.core,2)
    yhat=dotkron(yhat,ut{i})*reshape(TN.core{i},[TN.sz(i,1)*TN.sz(i,2),TN.sz(i,3)]); %O(dIr^2)
end

end

