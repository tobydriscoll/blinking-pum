function [ output ] = CoarseGlobalCorrection(rhs,tree)

cs = ASCoarseGlobalPreconditioner(tree,rhs);

z = rhs - LaplacianForward(tree,cs);

output = cs + ASPreconditioner(tree,z);
end

