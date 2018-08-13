function [ u,normres,normstep,numgm ] = PreconditionedNewton(f,Jac,init,PUApprox,tol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
normres = []; normstep = [];  numgm = [];

u = init;

% solve for the new value using plain Newton
for k = 1:100
    % evaluate the local corrections/solve local nonlinear problems
    [z,L,U,p] = ParPreconditionedNewtonForward(u,PUApprox,f,Jac);
    normres(k) = norm(z);
    
    if normres(k) < tol, break, end
    
    % find overall Newton step by GMRES
    tol_g = min(0.1,1e-10*norm(u)/normres(k));
    
    if 0 == tol_g
        tol_g = 0.1;
    end
    
    [s,~,~,~,gmhist] = gmres(@(x)JacobianFowardLU(PUApprox,L,U,p,x),-z,[],tol_g,100);
    
    normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
    
    % update
    u = u+s;
end

end

