function [ output ] = LaplacianForward(Tree,sol)

Tree.sample(sol);

LEAVES = Tree.collectLeaves();

output = zeros(size(sol));

step = zeros(length(LEAVES),1);

for k=2:length(LEAVES)
    step(k) = step(k-1) + length(LEAVES{k-1});
end

for k=1:length(LEAVES)
    
    dim = LEAVES{k}.degs;
    
    sol_k = sol(step(k)+(1:prod(dim)));
    
    lap = LEAVES{k}.linOp*sol_k;
    
    [~,~,in_border,~] = FindBoundaryIndex2DSides(dim,LEAVES{k}.domain,LEAVES{k}.outerbox);
    
    lap(in_border) = lap(in_border)-LEAVES{k}.Binterp*sol;
    
    output(step(k)+(1:length(LEAVES{k}))) = lap;
end

end


