function [un] = basisvectors(featurez,n,m)

%[un] = basisvectors(featurez,n,m)
% -------------
% Construct basis vectors from data samples
% 
% un        =   B-spline basis vectors
% 
% featurez  = Matrix of training data, rows are samples, columns are different features
% 
% n         = degree of B-spline
% 
% m         = number of knots intervals

[N, d]=size(featurez); 

In= m+n;
M = basismat(n);

knotdist = 1/m;         %distance between knots

indexes = floor(featurez/knotdist)+1;       %Calculate in which knot intervals the data samples fall.
indexes(indexes>m)= m;                      %Correction for when a data sample equals the upper limit (one).    
indexes(indexes<1)= 1;

inputs = (featurez/knotdist)-indexes+1;     %Map the datapoints to the knot intervals.

% The complexity of this operation is O(N*d*(n+1)^2) = O(N*d*n^2)
for i=1:d
bn = inputs(:,i).^[n:-1:0]*M;   % Construct the nonzero elements of the b-spline basis vectors using the matrix form.

un{i} = zeros(N,In);
for ii=1:N
   un{i}(ii,indexes(ii,i):indexes(ii,i)+n) = bn(ii,:); %Store them in the correct location within the basis vector. 
end
end

end

