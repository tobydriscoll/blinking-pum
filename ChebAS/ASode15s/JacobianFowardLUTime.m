% JacobianFowardLU
% Matrix free evaluation for the Jacobian of the SNK method.
%
% INPUT:
%    PUApproxArray: cell array of PUApprox approximations    
%            L,U,p: cell array of the LU decomposition of the local Jacobians. 
%                x: solution
%
% NOTE x is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
function [ output ] = JacobianFowardLUTime(PUApproxArray,L,U,p,x)

num_sols = length(PUApproxArray);

num_leaves = length(PUApproxArray{1}.leafArray);

sol_lengths = zeros(num_sols,1);

for i=1:num_sols
    sol_lengths(i) = length(PUApproxArray{i});
end

x = mat2cell(x,sol_lengths);

x_unpacked = x;



for i=1:num_sols
    
    if PUApproxArray{i}.is_packed
        
        %pull the boundry info for the packed functions. The boundary
        %info is stored within the tree itself.
        x_unpacked{i} = PUApproxArray{i}.Getunpackedvalues(x{i});
        
    else
        
        x_unpacked{i} = x{i};
    
    end
    
end

%x = reshape(x,length(PUApprox),num_sols);
[x,lens] = unpackPUvecs(cell2mat(x),PUApproxArray);

z = x;
%non parallel part
for k=1:num_leaves

    z_loc = mat2cell(zeros(size(x{k})),lens{k});
        
    for j=1:num_sols
        in_border = PUApproxArray{j}.leafArray{k}.inner_boundary;
        z_loc{j}(in_border) = PUApproxArray{j}.leafArray{k}.Binterp*x_unpacked{j};
    end
    
    z{k} = cell2mat(z_loc);
    
end

output = x;

%Parallel part
for k=1:num_leaves
    
    output{k} = U{k}\(L{k}\z{k}(p{k})) - x{k};
    
end

output = packPUvecs(output,PUApproxArray);