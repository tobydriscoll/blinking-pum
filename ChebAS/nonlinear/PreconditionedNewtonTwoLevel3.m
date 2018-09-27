function [ u,normres,normstep,numgm ] = PreconditionedNewtonTwoLevel3(f,Jac,init,PUApprox,Leaf,G,tol,tol_c,j)
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
    [z,L,U,p,Jac_hat] = ParPreconditionedTwoLevel3(u,PUApprox,Leaf,G,f,Jac,tol,tol_c);
    
    normres(k) = norm(z);
    
    normres(k)
    
    if normres(k) < tol(1), break, end
    
    
    % find overall Newton step by GMRES
    tol_g = min(0.1,norm(u)/norm(z)*1e-10);    
    
    [FJv,FJv_hat] = ComputeJac3(PUApprox,Leaf,G,Jac,u);
    
%    FJv_hat = blkdiag(FJv_hat{:});
    
    [L_hat,U_hat,p_hat] = lu(Jac_hat,'vector');
    J_hat.L = L_hat; J_hat.U = U_hat; J_hat.p = p_hat;
    
    [s,~,~,~,gmhist] = gmres(@(w)JacobianFoward2Level3(PUApprox,Leaf,G,L,U,p,J_hat,FJv,FJv_hat,w),-z,[],tol_g,80);
    
    normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
    
    numgm(k)
    
    % update
    u = u+s;
end

end

