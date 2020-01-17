function [TN,Vm,Vp] = initTT(un,r,d)

% [TN,Vm,Vp] = initTT(un,r,d)
% -------------
% Initializes a Tensor Train with random values.
%
% TN         =   initialized Tensor Train
% 
% Vm        =   cel containing left contracted wings of the TT
% 
% Vp        =   cel containing right contracted wings of the TT
% 
% d         =   number of dimensions/cores
%
% r         =   vector with TT-ranks
%
% un        =   B-spline basis vectors
% 

[N, ~]=size(un{1}); 

for i=1:d
    I(i)=size(un{i},2);     %size of the dimensions
end

% Randomly initialize TT-cores 
init=cell(1,d);
    for i=1:d
        init{i}.core=randn(r(i),I(i),r(i+1));
        init{i}.core=(init{i}.core/(norm(init{i}.core(:)))); %normalize cores to unit norm
        init{i}.r(1)=r(i);
        init{i}.I = I(i);
        init{i}.r(2)=r(i+1);
    end

% Store the core and its dimensions in a struct
for i=1:d
    TN.core{i}=init{i}.core;
    TN.sz(i,:)=[init{i}.r(1),init{i}.I,init{i}.r(2)];
end

Vp=cell(1,d);   
Vm=cell(1,d);   %cel containing left contracted wings of the TT

Vm{1}=ones(N,1);
Vp{d}=ones(N,1);

% initialize right-orthonormal cores with prescribed TN ranks
for i=d:-1:2
    Vp{i-1}=dotkron(Vp{i},un{i})*reshape(permute(TN.core{i},[3 2 1]),[r(i+1)*I(i),r(i)]); 
end


end

