function [ output ] = LaplacianForward(Tree,domain,sol)

Tree.sample(sol);

LEAVES = Tree.collectLeaves({});

output = [];

step_n = 0;

for k=1:length(LEAVES)
    
    dim = LEAVES{k}.degs;
    
    sol_k = sol(step_n+(1:prod(dim)));
    
    lap = LEAVES{k}.linOp*sol_k;
    
    
    
%     [out_border, in_border] = FindBoundaryGridIndex2D(dim,LEAVES{k}.domain(),domain);
%     grids = LEAVES{k}.leafGrids();
%     for j=1:4
%         if any(in_border{j,1}) && any(in_border{j,2})
%             approx = Tree.evalfZoneGrid({grids{1}(in_border{j,1}) grids{2}(in_border{j,2})});
%             [X_in,Y_in] = ndgrid(in_border{j,1},in_border{j,2});
%             IND = X_in & Y_in;
%             IND = IND(:);
%             lap(IND) = lap(IND)-approx(:);
%         end
%     end
    
    [out_border, in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
    lap(in_border) = lap(in_border)-LEAVES{k}.Binterp*sol;
    output = [output;lap];
    
    step_n = step_n + prod(dim);
    
end

end


