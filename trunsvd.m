function [Utrun,Strun,Vtrun,r] = trunsvd(C,e,Atest)

%[Utrun,Strun,Vtrun,r] = trunsvd(C,e,Atest)
% -------------
% Performs a truncated svd on matrix C for a given error e
% 
% C         = Matrix
% 
% e         = error
% 
% Atest     = Norm of original tensor

[U,S,V] = svd(C);
sings = diag(S).^2;
cumuls = cumsum(sings)/(Atest);
r = find(cumuls>1-e, 1);

Utrun = U(:,1:r);
Strun = S(1:r,1:r);
Vtrun = V(:,1:r);

end