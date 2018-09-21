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
function [z,l,u,p] = ParPreconditionedNewtonForwardTime(t,sol,rhs,PUApprox,evalF,hinvGak,Jac,M)

num_sols = length(sol)/length(PUApprox);

step = zeros(length(PUApprox.leafArray),1);

sol = reshape(sol,length(PUApprox),num_sols);
rhs = reshape(rhs,length(PUApprox),num_sols);

%Figure out starting index for each patch
for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end

for k=1:length(PUApprox.leafArray)
    
    degs = PUApprox.leafArray{k}.degs;
    
    %This will be sol_length*num_sols
    sol_loc{k} = sol(step(k)+(1:prod(degs)),:);
    rhs_loc{k} = rhs(step(k)+(1:prod(degs)),:);
    
    %This function returns the logical indicies of the gamma and outer
    %boundry interface. Out put is given for all indicies, as well as the
    %indicies along each of the sides
    [~,~,in_border{k},~] = FindBoundaryIndex2DSides(degs,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox);
    
    
    %This will be (interface length)*num_sols
    diff{k} = PUApprox.leafArray{k}.Binterp*sol;
    
end

%parallel step

leafs = PUApprox.leafArray;

for k=1:length(leafs)
    
    [z{k},l{k},u{k},p{k}] = local_inverse(leafs{k},sol_loc{k},t,rhs_loc{k},in_border{k},diff{k},evalF,hinvGak,num_sols,Jac,M{k});
    
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
function [c,l,u,p] = local_inverse(approx,sol_k,t,rhs_k,border_k,diff_k,evalF,hinvGak,num_sols,Jac,M)

%The residul is F(sol_k+z_k)
%            sol_k(border_k)+z_k(border_k)-B_k*u
%            (iterfacing at the zone interface)
    function [F] = residual(z)
        
        sol_length = length(approx);
        
        F = hinvGak*evalF(t,z+sol_k(:),approx)+rhs_k(:)-M*z;
        
        F = reshape(F,sol_length,num_sols);
        
        z = reshape(z,sol_length,num_sols)+sol_k;
        
        F(border_k,:) = z(border_k,:) - diff_k;
        
        F = F(:);
        
    end

    function J = jac_fun(z)
        
        sol_length = length(approx);
        
        J = hinvGak*Jac(t,z+sol_k(:),approx)-M;
        
        E = eye(sol_length);
        
        for i=1:num_sols
            ind = false(sol_length*num_sols,1);
            ind((i-1)*sol_length+(1:sol_length)) = border_k;
            
            J(ind,:) = zeros(sum(ind),num_sols*sol_length);
            J(ind,(i-1)*sol_length+(1:sol_length)) = E(border_k,:);
        end
    end


options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',1000,'FunctionTolerance',1e-4,'Display','off');
[c,~,~,~,~] = fsolve(@(u)sol_and_jac(@residual,@jac_fun,u),zeros(numel(sol_k),1),options);
c = c(:,end);

J = jac_fun(c);

[l,u,p] = lu(J,'vector');

%params = [20,-1,.5,0];
%tol = [1e-4 1e-4];
%[c,l,u,p] = nsoldAS(zeros(numel(sol_k),1),@residual,@jac_fun,tol,params);

%c = s(:,end);
end


