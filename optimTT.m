function [TN,Vm,Vp,res1,res2] = optimTT(TN,Vm,Vp,un,zeta,MAXITR,nselect,lambda,difforder)

% [TN,Vm,Vp,res1,res2] = optimTT(TN,Vm,Vp,un,zeta,MAXITR,nselect,lambda,gamma,difforder)
% -------------
% Optimizes the Tensor Train to find the B-spline weights.
%
% TN         =  Tensor Train representing the B-spline weights
%
% Vm         =  Cel containing right contracted wings of the TT
% 
% Vp         =  Cel containing left contracted wings of the TT
% 
% res1       =  Mean squared error on the training set
% 
% res2       =  Norm of the Tensor Train
% 
% un         =  B-spline basis vectors of the training set
% 
% zeta       = Outputs of the training set
% 
% MAXITR     = Maximum number of sweeps before termination
% 
% nselect    = Vector with batchsizes indicating the number of samples to use for updating the TT cores at each sweep
%           
% lambda     = Roughness penalty. Higher values penalizes the smoothness of adjacent weights, resulting in a smoother surface
% 
% gamma      = Variance relaxation, higher values allow more bias in the estimation of the weights. Currently unused.
% 
% difforder  = Order of difference to penalize. E.g, Second order penalizes the second derivative of the resulting curve.

res1=[];
res2=[];
itr=1;                          % counts number of iterations
ltr=1;                          % flag that checks whether we sweep left to right
sweepindex=1;                   % index that indicates which TT core will be updated

[N, ~]=size(un{1}); 
d=size(TN.sz,(1)); 
I = TN.sz(:,2);
r = [TN.sz(:,1);1];
P = diff(eye(I(1)),difforder);
PP = P'*P;


if isempty(nselect)
    nselect = ones([1,MAXITR])*N;
end

tic
while (itr <= MAXITR )
%     ---------------------updateTT-----------------;
        %Select random batch of data
        dataselect = randperm((N),nselect(itr));  
        
        %Construct regressor A
        A=dotkron(Vm{sweepindex}(dataselect,:),un{sweepindex}(dataselect,:),Vp{sweepindex}(dataselect,:));

        %Construct difference penalty matrices
        W = penalmat(TN,sweepindex,d,P,PP);
        
        %Sum the difference penalty matrices to penalize each dimension equally
        WWW = W{1};
        for s =2:d
            WWW = WWW + W{s};
        end
          
        %Solve linear subsystem 

             g=pinv(A'*A + (nselect(itr)/N)*lambda*WWW)*(A'*zeta(dataselect,:));

        %Update cores
        if ltr
%             left-to-right sweep, generate left orthogonal cores and update Vm
            [Q,R]=qr(reshape(g,[r(sweepindex)*(I(sweepindex)),r(sweepindex+1)])); 
            TN.core{sweepindex}=reshape(Q(:,1:r(sweepindex+1)),[r(sweepindex),I(sweepindex),r(sweepindex+1)]);
            TN.core{sweepindex+1}=reshape(R(1:r(sweepindex+1),:)*reshape(TN.core{sweepindex+1},[r(sweepindex+1),(I(sweepindex+1))*r(sweepindex+2)]),[r(sweepindex+1),I(sweepindex+1),r(sweepindex+2)]);
            Vm{sweepindex+1}=dotkron(Vm{sweepindex},un{sweepindex})*reshape(TN.core{sweepindex},[r(sweepindex)*I(sweepindex),r(sweepindex+1)]); 
        else
%             right-to-left sweep, generate right orthogonal cores and update Vp
            [Q,R]=qr(reshape(g,[r(sweepindex),(I(sweepindex))*r(sweepindex+1)])'); 
            TN.core{sweepindex}=reshape(Q(:,1:r(sweepindex))',[r(sweepindex),I(sweepindex),r(sweepindex+1)]);
            TN.core{sweepindex-1}=reshape(reshape(TN.core{sweepindex-1},[r(sweepindex-1)*(I(sweepindex-1)),r(sweepindex)])*R(1:r(sweepindex),:)',[r(sweepindex-1),I(sweepindex-1),r(sweepindex)]);
            Vp{sweepindex-1}=dotkron(Vp{sweepindex},un{sweepindex})*reshape(permute(TN.core{sweepindex},[3 2 1]),[r(sweepindex+1)*I(sweepindex),r(sweepindex)]);  
        end
    
        %Update sweep
     
                if ltr
                    sweepindex=sweepindex+1;
                    if sweepindex== d                
                        ltr=0;
                        timer=toc;
                    end
                else
                    sweepindex=sweepindex-1;
                    if sweepindex== 1                
                        ltr=1;
                        timer=toc;
                    end
                end
    
    
        % check residual after sweep
        if (sweepindex==d) || (sweepindex==1) % half a sweep
            
            res1(itr)=(1/nselect(itr))*norm(A*g-zeta(dataselect,:))^2; % check residual
            res2(itr)=(1/N)*lambda*(g'*WWW*g);

%             disp(["iteration:" itr timer])
%              disp(sweepindex)
            itr=itr+1; %update iteration
        end   
end  


end


function W = penalmat(TN,sweepindex,d,P,PP)

        for j=1:d
        Csize = TN.sz(j,:);
        
        Dm = reshape(permute(TN.core{j}, [2 1 3]), [Csize(2) Csize(1)*Csize(3)]);
        mDDm = reshape(Dm'*Dm, [Csize(1) Csize(3) Csize(1) Csize(3)]);                      %O(I*r^4)
        DD{j} = reshape(permute(mDDm,[1 3 2 4]), [Csize(1)*Csize(1) Csize(3)*Csize(3)]);
        PD = P*Dm;                                                                          %O(I^2*r^2)
        DPPD = reshape(PD'*PD, [Csize(1) Csize(3) Csize(1) Csize(3)]);                      %O(I^2*r^4)
        DWD{j} = reshape(permute(DPPD,[1 3 2 4]), [Csize(1)*Csize(1) Csize(3)*Csize(3)]);   
        eyez{j}= reshape(eye(Csize(1)), [Csize(1)^2 1]);                    
        eyep{j}= reshape(eye(Csize(2)), [Csize(2)^2 1]);                    
        end
        eyez{d+1}=1;
        
        for p= 1:d   %O(d^2*r^4)
            Dsize = TN.sz(sweepindex,:);   
                if sweepindex==p
                    D1 = eyez{sweepindex};
                    D2 = PP(:);
                    D3 = eyez{sweepindex+1};
                elseif sweepindex<p
                    D1 = eyez{sweepindex};  
                    D2 = eyep{sweepindex};   
                    D3= DWD{p}*eyez{p+1};           
                    for it=(p-1):-1:(sweepindex+1)               
                        D3 = DD{it}*D3;                     %O(d*r^4)
                    end                               
                elseif sweepindex>p   
                    D1= eyez{p}'*DWD{p};
       
                    for it=(p+1):(sweepindex-1)               
                        D1 = D1*DD{it};                     %O(d*r^4)
                    end
                    D1=D1';
                    D2 = eyep{sweepindex};   
                    D3= eyez{sweepindex+1}; 
                end
          
            WW = kron(D3, kron(D2, D1)); %O(I^4*r^4)
            Wtemp = permute(reshape(WW, [Dsize(1) Dsize(1) Dsize(2) Dsize(2) Dsize(3) Dsize(3)]), [1 3 5 2 4 6]);
            W{p} = reshape(Wtemp, [prod(Dsize) prod(Dsize)]); 
        end

end