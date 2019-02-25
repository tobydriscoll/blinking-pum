% ComputeJacs
% Computes Jacobian for NSK method
%
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
function [J,l,u,p] = ComputeJacs(sol,PUApproxArray,NonLinOps)

if ~iscell(PUApproxArray)
    PUApproxArray = {PUApproxArray};
end

%PUApprox.sample(sol);

num_sols = length(PUApproxArray);

num_leaves = length(PUApproxArray{1}.leafArray);

sol_lengths = zeros(num_sols,1);

for i=1:num_sols
    sol_lengths(i) = length(PUApproxArray{i});
end

%if any are packed, set boundary terms.
setBoundary(NonLinOps,PUApproxArray);
        

%Take [u1;u2;v1;v2] to {[u1;u2],[v1;v2]}
[ sol_loc,lens ] = unpackPUvecs(sol,PUApproxArray);

start_index = zeros(num_sols,1);

border = cell(num_leaves,num_sols);

for k=1:num_leaves

    for i=1:num_sols
        border{k}{i} = PUApproxArray{i}.leafArray{k}.inner_boundary;
    end
    
    start_index = start_index + lens{k};
    
end

%parallel step

for k=1:num_leaves
    
    [J{k},l{k},u{k},p{k}] = local_Jac(sol_loc{k},NonLinOps{k},border{k},lens{k});
    
end

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
    function [J,L,U,p] = local_Jac(sol_k,NonLinOps_k,border_k,lens_k)
        
        num_sols = length(lens_k);
        
        J = NonLinOps_k.jac(sol_k);

        index = 0;
        
        %This is supposed to account for the interfacing
        for i=1:num_sols
            
            E = eye(lens_k(i));
            
            total_length = sum(lens_k);
            
            ind = false(total_length,1);
            
            local_ind = index+(1:lens_k(i));
            
            ind(local_ind) = border_k{i};
            
            J(ind,:) = zeros(sum(ind),total_length);
            J(ind,local_ind) = E(border_k{i},:);
            
            index = index+lens_k(i);
        end
        
        [L,U,p] = lu(J,'vector');
    end
    
    function setBoundary(NonLinOps,PUApproxArray)
    for i=1:length(PUApproxArray{1}.leafArray)
        NonLinOps{i}.SetBoundaryVals();
    end
end





