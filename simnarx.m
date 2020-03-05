function [simoutput] = simnarx(TN,input,output,inlags,outlags,n,m)

%[un] = basisvectors(featurez,n,m)
% -------------
% Construct basis vectors from data samples
% 
% TN        =   Tensor Train struct
% 
% input  = input sequence
% 
% output  = output sequence
% 
% inlags = lags for the input
% 
% outlags = lags for the outlags
% 
% n         = degree of B-spline
% 
% m         = number of knots intervals

siminput = input;
simoutput = output;
beginz = max(outlags(end),inlags(end))+1;
for j = beginz:length(input)
simfeaturez = [simoutput(j-outlags(end:-1:2))' siminput(j-inlags(1:end))'];
simfeaturez(simfeaturez<0)=0+eps;
simfeaturez(simfeaturez>1)=1-eps;
simoutput(j) = evalspline(TN,simfeaturez,n,m);
end

end

