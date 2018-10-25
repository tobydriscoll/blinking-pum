function [ u,normres,normstep,numgm ] = NSKsolver(f,Jac,init,PUApprox,tol_n)
% NKSsolver
% The Newton Schwarz Krylov (NKS) solves nonlinear PDEs by
% preconditioning the linear solve for the Newton step. With method, the
% PDE is solved for overlaping subdomains with independent grids.
%
%
% INPUT:     
%             f: nonlinear residual function f(x,p) for solution x and local
%             approximation p. f(x,p) evaluates the residual on the domain
%             of p.
%
%            Jac: Jacobian function Jac(x,p) for solution x and local
%             approximation p. Jac(x,p) evaluates the residual on the domain
%             of p.
%
%           init: initial guess for all subdomains.
%
%       PUApprox: PUApprox approximation   
%
%          tol_n: [rel_tol, abs_tol] relative and absolute tolerance used
%                 for Newtons method.
%
% OUTPUT:
%              u: solution of PDE
%        normres: norm of SNK forward evaluation of each iteration
%       normstep: norm of Newton step for each iteration
%          numgm: number of GMRES iterations be iteration
%       normresf: norm of PDE residual of each iteration
%
% NOTE u,init is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
normres = []; normstep = [];  numgm = [];

u = init;


% solve for the new value using plain Newton
for k = 1:20
    
     % evaluate the local corrections/solve local nonlinear problems
    z = ParResidual(u,PUApprox,f);
    
    if k==1
        stop_tol = norm(z)*tol_n(1)+tol_n(2);
    end
    
    normres(k) = norm(z);
    
    normres(k)
    
    if normres(k) < stop_tol, break, end
    
    % find overall Newton step by GMRES
    tol_g = min(0.1,norm(u)/norm(z)*1e-10);
    
    if 0 == tol_g
        tol_g = 1e-10;
    end
    
    tol_g = 1e-10;
    
    [J,L,U,p] = ComputeJacs(u,PUApprox,Jac);  
    
    [s,~,~,~,gmhist] = gmres(@(x)ParLinearResidual(x,PUApprox,J),-z,[],tol_g,300,@(x)ASPreconditionerMultSols(PUApprox,U,L,p,x));
    
    normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
    
    numgm(k)
     
    % update
    u = u+s;
end

end

