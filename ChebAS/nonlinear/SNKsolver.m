function [ u,normres,normstep,numgm,normresf ] = SNKsolver(f,Jac,init,PUApprox,tol_n)
% SNKsolver
% The Schwarz Newton Krylov  method (SNK) solves nonlinear PDEs by nonlinear
% preconditioning Newton's method via a alternating Schwarz process.
% With this method, the PDE is solved for overlaping subdomains with 
% independent grids.
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
%                 for Newtons method, both with SNK and the local nonlinear
%                 solves.
%
% OUTPUT:
%              u: solution of PDE.
%
%        normres: norm of SNK forward evaluation of each iteration.
%
%       normstep: norm of Newton step for each iteration.
%
%          numgm: number of GMRES iterations be iteration.
%
%       normresf: norm of PDE residual of each iteration.
%
% NOTE u,init is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
normres = []; normstep = [];  numgm = []; normresf = [];

u = init;

% solve for the new value using plain Newton
for k = 1:20

    normresf(k) = norm(ParResidual(u,PUApprox,f));
    
    % evaluate the local corrections/solve local nonlinear problems
    [z,L,U,p] = SNK_forward_eval(u,PUApprox,f,Jac,tol_n);
    
    normres(k) = norm(z);
    
    
    if k==1
        stop_tol = norm(z)*tol_n(1)+tol_n(2);
    end
    
    normres(k)
    normresf(k)
    
    if normres(k) < stop_tol, break, end

   % tol_g = 1e-10;
   tol_g = min(max(1e-10,1e-10/normres(k)*norm(u)),1e-1);
   tol_g 
    %solve the newton step
    [s,~,~,~,gmhist] = gmres(@(x)JacobianFowardLU(PUApprox,L,U,p,x),-z,[],tol_g,200);
    
    normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
    
    numgm(k)
    
    % update
    u = u+s;
end

end

