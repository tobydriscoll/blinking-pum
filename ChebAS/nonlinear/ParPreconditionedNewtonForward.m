% INPUT:
%      PUApprox: PUApprox approximation
%         sol: given solution
%       evalF: residual function which returns Jacobian
%    num_sols: number of solutions
%
% OUTPUT:
%          z: correction of solution
%          J: cell array of local Jacobians
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [z,l,u,p,J] = ParPreconditionedNewtonForward(sol,PUApprox,evalF,Jac,tol_n)

num_sols = length(sol)/length(PUApprox);

step = zeros(length(PUApprox.leafArray),1);

sol = reshape(sol,length(PUApprox),num_sols);

%Figure out starting index for each patch
for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end

for k=1:length(PUApprox.leafArray)
    
    degs = PUApprox.leafArray{k}.degs;
    
    %This will be sol_length*num_sols
    sol_loc{k} = sol(step(k)+(1:prod(degs)),:);

    
    %This function returns the logical indicies of the gamma and outer
    %boundry interface. Out put is given for all indicies, as well as the
    %indicies along each of the sides
    [~,~,in_border{k},~] = FindBoundaryIndex2DSides(degs,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox);
    
    
    if ~PUApprox.iscoarse
        diff{k} = PUApprox.leafArray{k}.Binterp*sol;
    else
        diff{k} = PUApprox.leafArray{k}.CBinterp*sol;
    end
    
end

%parallel step

leafs = PUApprox.leafArray;

for k=1:length(leafs)
    
    [z{k},l{k},u{k},p{k},J{k}] = local_inverse(leafs{k},sol_loc{k},in_border{k},diff{k},evalF,num_sols,Jac,tol_n);
    z{k} = reshape(z{k},length(leafs{k}),num_sols);
end

z = cell2mat(z');
z = z(:);

end

% INPUT:
%      approx: leaf approximation
%       sol_k: given solution
%    border_k: border index for interface
%      diff_k: precomputed interface zone interpolation
%       evalF: residual function which returns Jacobian
%    num_sols: number of solutions
%  sol_length: length of solution
%
% OUTPUT
%           c: correction of solution
%          Jk: local Jocabian
function [c,l,u,p,J] = local_inverse(approx,sol_k,border_k,diff_k,evalF,num_sols,Jac,tol_n)

%The residul is F(sol_k+z_k)
%            sol_k(border_k)+z_k(border_k)-B_k*u
%            (iterfacing at the zone interface)
    function [F] = residual(z)
        
        sol_length = length(approx);
        
        F = evalF(sol_k(:)-z,approx);
        
        F = reshape(F,sol_length,num_sols);
        
        z = reshape(z,sol_length,num_sols);
        
        F(border_k,:) = sol_k(border_k,:)-z(border_k,:) - diff_k;
        
        F = F(:);
        
    end

    function J = jac_fun(z)
        
        sol_length = length(approx);
        
        J = -Jac(sol_k(:)-z(:),approx);
        
        E = eye(sol_length);
        
        for i=1:num_sols
            ind = false(sol_length*num_sols,1);
            ind((i-1)*sol_length+(1:sol_length)) = border_k;
            
            J(ind,:) = zeros(sum(ind),num_sols*sol_length);
            J(ind,(i-1)*sol_length+(1:sol_length)) = -E(border_k,:);
        end
        
    end

tol_c = tol_n(1)*norm(residual(zeros(numel(sol_k(:)),1)))+tol_n(2);
options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',500,'TolFun',tol_c,'Display','iter');

c = fsolve(@(u)sol_and_jac(@residual,@jac_fun,u),zeros(numel(sol_k),1),options);

c = c(:,end);

%[out_border_s,~,~,~] = FindBoundaryIndex2DSides(approx.degs,approx.domain,approx.outerbox);

%north = out_border_s{4};

%params = [200,-1,.5,0];
%tol = [1e-3 1e-2];

%c2 = nsoldAS(sol_k(:),@residual,@jac_fun,tol_n,params);

%norm_c = norm(residual(c));

%R = residual(c);

J = -jac_fun(c);

%AJ = jacobi(@residual,c);

[l,u,p] = lu(J,'vector');

end




