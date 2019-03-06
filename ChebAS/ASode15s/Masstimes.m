function R = Masstimes(PUApproxArray,M,y)

[y_loc,lens] = unpackPUvecs(y,PUApproxArray);

num_leaves = length(PUApproxArray{1}.leafArray);

num_sols = length(PUApproxArray);

for i=1:num_leaves
    
    R_i = M{i}*y_loc{i};
    R_i = mat2cell(R_i,lens{i});
    
    for j=1:num_sols
       R_i{j}(PUApproxArray{j}.leafArray{i}.inner_boundary) = 0;
    end
    
    R{i} = cell2mat(R_i);
    
end

R = packPUvecs(R,PUApproxArray);

end

