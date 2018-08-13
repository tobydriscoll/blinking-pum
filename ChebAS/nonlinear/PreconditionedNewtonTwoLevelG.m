function [ u,normres,normstep,numgm ] = PreconditionedNewtonTwoLevelG(f,Jac,init,PUApprox,tol,j)
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
    [z,L,U,p,g_mat,g_p_mat,J_er] = ParPreconditionedTwoLevelG(u,PUApprox,f,Jac,j);
    
    normres(k) = norm(z);
    
    normres(k)
    
    if normres(k) < tol, break, end
    
    % find overall Newton step by GMRES
    tol_g = min(0.1,1e-10*norm(u)/normres(k));    
    
    [s,~,~,~,gmhist] = gmres(@(w)JacobianFoward2LevelG(PUApprox,L,U,p,g_mat,g_p_mat,J_er,T_hat,w,j),-z,[],tol_g,100);
    
    normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
    
    % update
    u = u+s;
end

end

