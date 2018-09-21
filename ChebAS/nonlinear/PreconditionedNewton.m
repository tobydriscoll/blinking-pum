function [ u,normres,normstep,numgm ] = PreconditionedNewton(f,Jac,init,PUApprox,tol,tol_n,Leaf,G)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
normres = []; normstep = [];  numgm = [];

u = init;

% solve for the new value using plain Newton
for k = 1:100

    normres(k) = norm(ParResidual(u,PUApprox,f));
    
    normres(k)

    if normres(k) < tol, break, end
    
    % evaluate the local corrections/solve local nonlinear problems
    [z,L,U,p,J] = ParPreconditionedNewtonForward(u,PUApprox,f,Jac,tol_n);
    
    
    
    % find overall Newton step by GMRES
    tol_g = min(0.1,norm(u)/norm(z)*1e-10);
    
    if 0 == tol_g
        tol_g = 1e-10;
    end
    
    %[J,L,U,p] = ComputeJacs(u,PUApprox,Jac);
    
    fz = ParLinearResidualFun(z,PUApprox,J);
    
    
    Len = length(PUApprox);
    
    sol_hat = [];
    
%     num_sols = length(u)/Len;
%     
%     for i=1:num_sols
%         PUApprox.sample( u((1:Len) + (i-1)*Len) - z((1:Len) + (i-1)*Len));
%         F = PUApprox.evalfGrid(G);
%         sol_hat = [sol_hat;F(:)];
%     end
    
    
    %Jac_hat = Jac(sol_hat,Leaf);
    
    
    [s,~,~,~,gmhist] = gmres(@(x)ParLinearResidual(x,PUApprox,J),-fz,[],tol_g,300,@(x)ASPreconditionerMultSols(PUApprox,U,L,p,x));
    
    %[s,~,~,~,gmhist] = gmres(@(x)ParLinearResidual(x,PUApprox,J),-fz,[],tol_g,300,@(x)CoarseGlobalCorrection(x,J,L,U,p,Jac_hat,PUApprox,Leaf,G));
    
    %[s,~,~,~,gmhist] = gmres(@(x)JacobianFowardLU(PUApprox,L,U,p,x),-z,[],tol_g,300);
    
    normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
    
    numgm(k)
    
    % update
    u = u+s;
end

end

