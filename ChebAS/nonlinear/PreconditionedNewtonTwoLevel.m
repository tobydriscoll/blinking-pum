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
    
    normres(k)
    
    if normres(k) < newn_tol, break, end
    
    % find overall Newton step by GMRES
    tol_g = min(0.1,tol(1)*norm(u)+tol(2));    
    
    [FJv,FJv_hat] = ComputeJac(PUApprox,Jac,u);
    
    FJv_hat = blkdiag(FJv_hat{:});
    
    [s,~,~,~,gmhist] = gmres(@(w)JacobianFoward2Level(PUApprox,L,U,p,J_v_pls_er,T_hat,FJv,FJv_hat,w,j),-z,[],tol_g,80);
    
    normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
    
    numgm(k)
    
    % update
    u = u+s;
end

end

