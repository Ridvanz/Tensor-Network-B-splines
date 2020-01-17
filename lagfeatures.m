function [featurez,zeta,tfeaturez,yt] = lagfeatures(input,tinput,output,toutput,inlags,outlags)

ending = length(input)-max(outlags(end),inlags(end))-1;
endingt = length(tinput)-max(outlags(end),inlags(end))-1;

% generate features using input lags
for l = 1:length(inlags)
u(:,l) = input(end-inlags(l)-ending:end-inlags(l));
uv(:,l) = tinput(end-inlags(l)-endingt:end-inlags(l));
end
% generate features using output lags
for l = 1:length(outlags)
y(:,l) = output(end-outlags(l)-ending:end-outlags(l));
yv(:,l) = toutput(end-outlags(l)-endingt:end-outlags(l));
end

% Features are lagged inputs and outputs
featurez = [y(:,end:-1:2)  u];
% Target is output with zero lag
zeta = y(:,1);

tfeaturez = [yv(:,end:-1:2)  uv];
yt = yv(:,1);

end

