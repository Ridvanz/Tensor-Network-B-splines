function [TN] = trunttsvd(A,rk)

% Tensor Train decomposition with svd, given a truncation rank r

Sizes = size(A);
d = ndims(A);
G  = cell(d,1);

C=A;

for i = 1:d-1
    
C = reshape(C,rk(i)* Sizes(i),[]);

[U,S,V] = svd(C);

Utrun = U(:,1:rk(i+1));
Strun = S(1:rk(i+1),1:rk(i+1));
Vtrun = V(:,1:rk(i+1));


G{i}= reshape(Utrun, [rk(i), Sizes(i), rk(i+1)]);
C = Strun*Vtrun';

end

G{d} = C;

for i=1:d-1
    TN.core{i}=G{i};
    TN.sz(i,:)=size(G{i});
end
TN.core{d}=G{d};
TN.sz(d,:)=[size(G{d}) 1];

end