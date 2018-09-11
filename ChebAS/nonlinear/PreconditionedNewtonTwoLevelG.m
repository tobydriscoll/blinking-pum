function [ u,normres,normstep,numgm ] = PreconditionedNewtonTwoLevelG(f,Jac,init,PUApprox,tol,tol_c,j)
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
    [z,L,U,p,g_mat,g_p_mat,J_er] = ParPreconditionedTwoLevelG(u,PUApprox,f,Jac,tol,tol_c,j);
    
    normres(k) = norm(z);
    
    normres(k)
    
    if normres(k) < tol, break, end
    
    % find overall Newton step by GMRES
    tol_g = min(0.1,tol(1)*norm(z)+tol(2));    
    
    [s,~,~,~,gmhist] = gmres(@(w)JacobianFoward2LevelG(PUApprox,L,U,p,g_mat,g_p_mat,J_er,T_hat,w,j),-z,[],tol_g,100);
    
    normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
    
    numgm(k)
    % update
    u = u+s;
end

end

