function R = InterFaceError(PUApproxArray,y)

[y_loc,lens] = unpackPUvecs(y,PUApproxArray);

num_leaves = length(PUApproxArray{1}.leafArray);

num_sols = length(PUApproxArray);

for i=1:num_leaves
    
    R_i = zeros(size(y_loc{i}));
    
    R_i = mat2cell(R_i,lens{i});
    
    y_i = mat2cell(y_loc{i},lens{i});
    
    for j=1:num_sols
       R_i{j}(PUApproxArray{j}.leafArray{i}.inner_boundary) = y_i{j}(PUApproxArray{j}.leafArray{i}.inner_boundary);
    end
    
    R{i} = cell2mat(R_i);
    
end

R = packPUvecs(R,PUApproxArray);

R = norm(R,inf);

end
