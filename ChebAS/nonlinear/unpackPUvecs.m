function [ sol_loc ] = unpackPUvecs(sol,PUApproxArray)

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

start_index = zeros(num_sols,1);

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

