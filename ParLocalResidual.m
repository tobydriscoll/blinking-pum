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
function [z] = ParLocalResidual(t,sol,dt,PUApproxArray,NonLinOps,alpha,set_interface)

if nargin<7
	set_interface = false;
end

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
		sol_unpacked{i} = PUApproxArray{i}.Getunpackedvalues(sol{i});
	else
		sol_unpacked{i} = sol{i};
	end
	
end

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

%parallel step
for k=1:num_leaves
	z{k} = local_resid(sol_loc{k},t,dt,diff{k},border{k},NonLinOps{k},lens{k},set_interface,alpha);
end

z = packPUvecs(z,PUApproxArray);

end

function F = local_resid(sol_k,t,dt,diff_k,border_k,NonLinOps_k,lens_k,set_interface,alpha)

num_sols = length(lens_k);
F = dt*NonLinOps_k.timederiv(t,sol_k);
F = mat2cell(F,lens_k);
sol_k_c =  mat2cell(sol_k,lens_k);
for i=1:num_sols
	F{i}(border_k{i}) = alpha*(sol_k_c{i}(border_k{i}) - diff_k{i});
end
F = cell2mat(F);

end