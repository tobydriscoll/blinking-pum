% INPUT:     
%      sol: given solution at patches
%      sol2: second given solution 
%      PUApprox: PUApprox approximation   
%      jac_f(x,leaf): function that returns jacobian given local solution x and
%      patch leaf.
%
% OUTPUT:
%          Mat: sparse matrix used for Jacobian of Coarse correction.
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [ J,Jc ] = ComputeJac(PUApprox,jac_f,sol)
%assume sol is the correct coarse length
num_sols = length(sol)/length(PUApprox);

step = zeros(length(PUApprox.leafArray),1);

sol = reshape(sol,length(PUApprox),num_sols);

%Figure out starting index for each patch
for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end

ii = [];
jj = [];
zz = [];

for k=1:length(PUApprox.leafArray)
    
    degs = PUApprox.leafArray{k}.degs;
    
    sol_loc = sol(step(k)+(1:prod(degs)),:);
     
    J{k} = ComptuteLocalJac(sol_loc,PUApprox.leafArray{k},jac_f,num_sols);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sol_loc_hat = [];
    
    for i=1:num_sols
        sol_loc_hat = [sol_loc_hat PUApprox.leafArray{k}.Fine2Coarse(sol_loc(:,i))];
    end

    PUApprox.leafArray{k}.Coarsen();
    
    Jc{k} = ComptuteLocalJac(sol_loc_hat,PUApprox.leafArray{k},jac_f,num_sols);
    
    PUApprox.leafArray{k}.Refine();

end

end

function J = ComptuteLocalJac(sol_loc,leaf,jac_f,num_sols)

    degs = leaf.degs;
    
    [~,~,in_border,~]  = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);

    J = jac_f(sol_loc(:),leaf);
        
    sol_length = prod(degs);
    
    E = eye(sol_length);
    
    for i=1:num_sols
        ind = false(sol_length*num_sols,1);
        ind((i-1)*sol_length+(1:sol_length)) = in_border;
        
        J(ind,:) = zeros(sum(ind),num_sols*sol_length);
        J(ind,(i-1)*sol_length+(1:sol_length)) = E(in_border,:);
    end
end





