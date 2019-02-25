function [ u,normres,normstep,numgm,tol_g ] = SNKsolver(init,PUApproxArray,NonLinOps,tol_n)
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
normres = []; normstep = [];  numgm = []; linres = [];

u = init;

% solve for the new value using plain Newton
for k = 1:20
    
    % evaluate the local corrections/solve local nonlinear problems
    [z,L,U,p] = SNK_forward_eval(u,PUApproxArray,NonLinOps);
    
    normres(k) = norm(z);
    
    
    if k==1
        stop_tol = norm(z)*tol_n(1)+tol_n(2);
    end
    
    normres(k)
    
    if normres(k) < stop_tol, break, end

    if k==1
        tol_g(k) = 1e-4;
    else
        %tol_g(k) = min(max(abs(normres(k)-linres(k-1))/normres(k-1),tol_g(k-1)^((1+sqrt(5))/2)),1e-2);
        tol_g(k) = max(min(tol_g(k-1),1e-4*(normres(k)/normres(k-1))^2),1e-10);
    end
    
    tol_g(k)
    
    [s,~,~,~,gmhist] = gmres(@(x)JacobianFowardLUTime(PUApproxArray,L,U,p,x),-z,[],tol_g(k),200);
    
    linres(k) = gmhist(end);
    
    normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
    
    numgm(k)
    
    % update
    u = u+s;
end

end

