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
function [J,l,u,p] = ComputeJacsTime(t,sol,PUApproxArray,NonLinOps,hinvGak,M,loc_sub_ind)

if ~iscell(PUApproxArray)
    PUApproxArray = {PUApproxArray};
end

take_sub_ind = nargin>6;

%PUApprox.sample(sol);

num_sols = length(PUApproxArray);

num_leaves = length(PUApproxArray{1}.leafArray);

sol_lengths = zeros(num_sols,1);

if ~take_sub_ind
   loc_sub_ind = cell(1, num_leaves);
end

for i=1:num_sols
    sol_lengths(i) = length(PUApproxArray{i});
end

%if any are packed, set boundary terms.
%setBoundary(NonLinOps,PUApproxArray);
        

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

        if iscell(M)
            [J{k},l{k},u{k},p{k}] = local_Jac(t,sol_loc{k},NonLinOps{k},border{k},lens{k},hinvGak,M{k},loc_sub_ind{k});
        else
            [J{k},l{k},u{k},p{k}] = local_Jac(t,sol_loc{k},NonLinOps{k},border{k},lens{k},hinvGak,M,loc_sub_ind{k});
        end
    
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
    function [J,L,U,p] = local_Jac(t,sol_k,NonLinOps_k,border_k,lens_k,hinvGak,M_k,loc_sub_ind)
        
        num_sols = length(lens_k);
        
        
        J = hinvGak*NonLinOps_k.jac(t,sol_k)-M_k;
        
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
        
        J2 = J;
        
        if ~isempty(loc_sub_ind)
            J2 = J2(loc_sub_ind,loc_sub_ind);
        end
        
        [L,U,p] = lu(J2,'vector');
    end
    
function setBoundary(NonLinOps,PUApproxArray)
    for i=1:length(PUApproxArray{1}.leafArray)
        NonLinOps{i}.SetBoundaryVals();
    end
end





