% INPUT:
%                sol_k: cell Array of local solutions 
%
%       PUApproxArrray: Cell Array of PUApprox approximations
%
% OUTPUT:
%
%                 sols: single vector of solutions.
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sols = [u1;u2;v1;v2].
function [ sols ] = packPUvecs(sol_k,PUApproxArray)

num_sols = length(PUApproxArray);

num_leaves = length(PUApproxArray{1}.leafArray);

sols = cell(num_sols,1);

for k=1:num_leaves
    
    index = 0;
    
    for i=1:num_sols
        
        loc_ind = index+(1:length(PUApproxArray{i}.leafArray{k}));
        sols{i} = [sols{i};sol_k{k}(loc_ind)];
        index = index + length(PUApproxArray{i}.leafArray{k});
    end
    
end

sols = cell2mat(sols);



