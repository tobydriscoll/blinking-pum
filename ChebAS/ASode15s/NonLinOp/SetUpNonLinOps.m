function [ NonLinOps ] = SetUpNonLinOps(PUApproxArray,rfun,jfun,boundf)

if ~iscell(PUApproxArray)
    PUApproxArray = {PUApproxArray};
end

num_sols = length(PUApproxArray);

num_leaves = length(PUApproxArray{1}.leafArray);

    for i=1:num_leaves
        
        for j=1:num_sols
            Chebs{j} = PUApproxArray{j}.leafArray{i};
        end
        
        NonLinOps{i} = NonLinOpIVP2D(Chebs,rfun,jfun,boundf);
        
    end
    
end

