function Mn = basismat(n)
%UNTITLED3 Construct the basis matrix for uniform B-splines of degree n 
%   Code is an implementation of the recursive method of the following
%   paper: General matrix representations for B-splines, - Kaihuai Qin

M{1} = 1;

for k = 2:(n+1)
    
D1 = diag(1:(k-1),0);
D2 = diag((k-2):-1:0,1);
D = [D1 zeros(k-1,1)] + D2(1:k-1,:);
M{k} = 1/(k-1) *( [M{k-1} ; zeros(1,k-1)] * D + [zeros(1,k-1) ;M{k-1}] * diff(eye(k)));
end

Mn = flipud(M{end});

end

