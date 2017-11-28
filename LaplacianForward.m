function [ output ] = LaplacianForward(Tree,domain,sol)

Tree.sample(sol);

LEAVES = Tree.collectLeaves({});

output = [];

step_n = 0;

for k=1:length(LEAVES)
    
    dim = LEAVES{k}.degs;
    
    sol_k = sol(step_n+(1:prod(dim)));
    
    lap = LEAVES{k}.linOp*sol_k;
    
    [out_border_c,out_border,in_border,in_border_c,in_border_g] = FindBoundaryIndex2DSides(dim,LEAVES{k}.domain(),LEAVES{k}.outerbox);
    
    for i=1:4
        if any(in_border_c{i})
            lap(in_border_c{i}) = lap(in_border_c{i})-LEAVES{k}.Binterp{i}*sol;
        end
    end
    
%    lap(in_border) = lap(in_border)-LEAVES{k}.Binterp*sol;
    
    output = [output;lap];
    
    step_n = step_n + prod(dim);
    
end

end


