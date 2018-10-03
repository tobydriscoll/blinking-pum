function [ u,normres,normstep,numgm ] = PreconditionedNewtonTwoLevel(f,Jac,init,PUApprox,newn_tol,tol,tol_c,j)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
normres = []; normstep = [];  numgm = [];

u = init;

num_sols = length(u)/length(PUApprox);

T_hat  = CoarseInterfaceInterp(PUApprox,num_sols);

% solve for the new value using plain Newton
for k = 1:100
        
    % evaluate the local corrections/solve local nonlinear problems
    [z,L,U,p,J_v_pls_er] = ParPreconditionedTwoLevel(u,PUApprox,f,Jac,j,tol,tol_c);
    
    normres(k) = norm(z);
    
        if k==1
            stop_tol = norm(z)* newn_tol(1) + newn_tol(2);
        end
    
    normres(k)
    
    if normres(k) < stop_tol, break, end
    
    
    [FJv,FJv_hat] = ComputeJac(PUApprox,Jac,u);
    
    FJv_hat = blkdiag(FJv_hat{:});
    
    [s,~,~,~,gmhist] = gmres(@(w)JacobianFoward2Level(PUApprox,L,U,p,J_v_pls_er,T_hat,FJv,FJv_hat,w,j),-z,[],1e-10,180);
    
    normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
    
    numgm(k)
    
    % update
    u = u+s;
end

end

