% INPUT:
%                  sol: vector of all solutions
%
%       PUApproxArrray: Cell Array of PUApprox approximations
%
% OUTPUT:
%
%              sol_loc: cell array of local solutions
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [ sol_loc,lens ] = unpackPUvecs(sol,PUApproxArray)

num_sols = length(PUApproxArray);

num_leaves = length(PUApproxArray{1}.leafArray);

sol_lengths = zeros(1,num_sols);

for i=1:num_sols
    sol_lengths(i) = length(PUApproxArray{i});
end

sol = mat2cell(sol,sol_lengths);

start_index = zeros(1,num_sols);

for k=1:num_leaves
    
    sol_loc{k} = [];
    lens{k} = [];
    PUApproxArray{i}.leafArray{k};
    
    
    for i=1:num_sols
        len = length(PUApproxArray{i}.leafArray{k});
        sol_loc{k} = [sol_loc{k};sol{i}(start_index(i)+(1:len))];
        
        lens{k} = [lens{k} len];
        
    end
    
    
    start_index = start_index + lens{k};
    

end

