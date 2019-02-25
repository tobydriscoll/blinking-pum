% ParResidual
% This method computes the residual of a PDE, as well as the boundary
% differences along the interface boundaries.
%
% INPUT:     
%           sol: given solution at patches
%
%      PUApprox: PUApprox approximation
%
%             f: nonlinear residual function f(x,p) for solution x and local
%             approximation p. f(x,p) evaluates the residual on the domain
%             of p.
%
% OUTPUT:
%          z: residual of solution, identity at inner boundary of patches
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [z] = ParResidual(sol,PUApproxArray,NonLinOps)

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

sol = mat2cell(sol,sol_lengths);

sol_unpacked = sol;

setBoundary(NonLinOps,PUApproxArray);
        
for i=1:num_sols
    
    if PUApproxArray{i}.is_packed
        

        %pull the boundry info for the packed functions. The boundary
        %info is stored within the tree itself.
        sol_unpacked{i} = PUApproxArray{i}.Getunpackedvalues(sol{i});
        
    else
        
        sol_unpacked{i} = sol{i};
    
    end
    
end

%Take [u1;u2;v1;v2] to {[u1;u2],[v1;v2]}
[ sol_loc,lens ] = unpackPUvecs(cell2mat(sol),PUApproxArray);

start_index = zeros(num_sols,1);

diff = cell(num_leaves,num_sols);
border = cell(num_leaves,num_sols);

for k=1:num_leaves

    for i=1:num_sols
        %This will be (interface length)*num_sols
        diff{k}{i} = PUApproxArray{i}.leafArray{k}.Binterp*sol_unpacked{i};
        border{k}{i} = PUApproxArray{i}.leafArray{k}.inner_boundary;
    end
    
    start_index = start_index + lens{k};
    
end

for k=1:num_leaves
        z_loc{k} = local_residual(sol_loc{k},diff{k},border{k},NonLinOps{k},lens{k});
end

%Take {[u1;u2],[v1;v2]} to [u1;u2;v1;v2]
z = packPUvecs(z_loc,PUApproxArray);

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
    function [F] = local_residual(sol_k,diff_k,border_k,NonLinOps_k,lens_k)
        
        num_sols = length(lens_k);
        
        F = NonLinOps_k.residual(sol_k);
        
        F = mat2cell(F,lens_k);
        
        sol_k_c =  mat2cell(sol_k,lens_k);
        
        for i=1:num_sols
            F{i}(border_k{i}) = sol_k_c{i}(border_k{i}) - diff_k{i};
        end
        
       F = cell2mat(F);
        
    end
    
function setBoundary(NonLinOps,PUApproxArray)
    for i=1:length(PUApproxArray{1}.leafArray)
        NonLinOps{i}.SetBoundaryVals();
    end
end

