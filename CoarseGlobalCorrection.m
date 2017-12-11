function [ output ] = CoarseGlobalCorrection(rhs,tree,domain)

cs = ASCoarseGlobalPreconditioner(tree,domain,rhs);

z = rhs - LaplacianForward(tree,domain,cs);

output = cs + ASPreconditioner(tree,domain,z);

LEAVES = tree.collectLeaves({});

% g = [];
% 
% tree.sample(output);
% 
% for i=1:length(LEAVES)
%     grid = LEAVES{i}.leafGrids();
%     G = tree.evalf(grid);
%     g = [g;G(:)];
% end
% 
% output = g;
end

