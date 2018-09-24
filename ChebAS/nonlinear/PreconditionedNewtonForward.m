function [ u,normres,normstep,numgm ] = PreconditionedNewtonForward(f,Jac,init,PUApprox,tol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
normres = []; normstep = [];  numgm = [];

u = init;

% solve for the new value using plain Newton
for k = 1:100
    
     % evaluate the local corrections/solve local nonlinear problems
    z = ParResidual(u,PUApprox,f);
    
    normres(k) = norm(z);
    
    normres(k)
    
    if normres(k) < tol(1), break, end
    
    
    % find overall Newton step by GMRES
    % find overall Newton step by GMRES
    tol_g = min(0.1,norm(u)/norm(z)*1e-10);
    
    if 0 == tol_g
        tol_g = 1e-10;
    end
        
    [J,L,U,p] = ComputeJacs(u,PUApprox,Jac);  
    
    [s,~,~,~,gmhist] = gmres(@(x)ParLinearResidual(x,PUApprox,J),-z,[],tol_g,300,@(x)ASPreconditionerMultSols(PUApprox,U,L,p,x));
    
    normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
    
    numgm(k)
     
    % update
    u = u+s;
end

end

