% ParLinearResidual
%
% This method computes the residual the linearized PDE
%
% INPUT:     
%      sol: given solution at patches
%      PUApprox: PUApprox approximation   
%      J: Jcell array of the local Jacobians.
%
% OUTPUT:
%          z: residual of solution, identity at inner boundary of patches
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [z] = ParLinearResidual(sol,PUApprox,J)

%PUApprox.sample(sol);

num_sols = length(PUApproxArray);

num_leaves = length(PUApproxArray{1}.leafArray);

sol_lengths = zeros(num_sols,1);

for i=1:num_sols
    sol_lengths(i) = length(PUApproxArray{i});
end

sol = mat2cell(sol,sol_lengths);

sol_unpacked = sol;

for i=1:num_sols
    
    if PUApproxArray{i}.is_packed
        
        %pull the boundry info for the packed functions. The boundary
        %info is stored within the tree itself.
        sol_unpacked{i} = PUApproxArray{i}.Getunpackedvalues(sol{i},true);
        
    else
        
        sol_unpacked{i} = sol{i};
    
    end
    
end

[ sol_loc,lens ] = unpackPUvecs(cell2mat(sol),PUApproxArray);

start_index = zeros(num_sols,1);

diff = cell(num_leaves,num_sols);
border = cell(num_leaves,num_sols);

%parallel step
for k=1:num_leaves

    for i=1:num_sols
        %This will be (interface length)*num_sols
        diff{k}{i} = PUApproxArray{i}.leafArray{k}.Binterp*sol_unpacked{i};
        border{k}{i} = PUApproxArray{i}.leafArray{k}.inner_boundary;
    end
    
    start_index = start_index + lens{k};
    
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
function F = local_residual(approx,sol_k,border_k,diff_k,J,num_sols)
     
        F = J*sol_k(:);  
        
        sol_length = length(approx);

        F = reshape(F,sol_length,num_sols);
        
        F(border_k,:) = sol_k(border_k,:) - diff_k;
              
end

