function [ u,normres,normstep,numgm ] = SNK2levelsolver(f,Jac,init,PUApprox,newn_tol,tol_c,j)
% SNKsolver
% The Schwarz Newton Krylov 2 level method (SNK2level) solves nonlinear 
% PDEs by nonlinear preconditioning the PDE with the a 2 level FAS method 
% using the alternating Schwarz process. With this method, the PDE is 
% solved for overlaping subdomains with independent grids.
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
%          tol_c: tolerance used for nonlinear coarse correction solve.
%
%              j: difference by power of two between coarse and fine grids.
%                 Grids are chosen from 2,5,9,33... . For example if we
%                 have coarse and fine grids of 9,33 (in each dimension)
%                 then j=2. (This should just be apart of PUApprox).
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
% NOTE u,init is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
normres = []; normstep = [];  numgm = [];

u = init;

num_sols = length(u)/length(PUApprox);

% solve for the new value using plain Newton
for k = 1:100
    
    % evaluate the local corrections/solve local nonlinear problems
    [z,L,U,p,J_v_pls_er,c_sol] = SNK2level_forward_eval(u,PUApprox,f,Jac,j,newn_tol,tol_c);
    
    normres(k) = norm(z);
    
    if k==1
        stop_tol = norm(z)* newn_tol(1) + newn_tol(2);
    end
    
    normres(k)
    
    if normres(k) < stop_tol, break, end
    
    
    [FJv,FJv_hat] = ComputeJac(PUApprox,Jac,u);
    
    [Lc,Uc,Pc,Qc] = lu(J_v_pls_er);

    tol_g = min(max(1e-10,1e-4/normres(k)*norm(u)),1e-1);
    
    tol_g = 1e-4;
    
    [s,~,~,~,gmhist] = gmres(@(w)JacobianFoward2Level(PUApprox,L,U,p,Lc,Uc,Pc,Qc,FJv,FJv_hat,w,j),-z,[],tol_g,600);
    
    normstep(k) = norm(s);  numgm(k) = length(gmhist) - 1;
    
    numgm(k)
    
    % update
    u = u+s;
end

end

