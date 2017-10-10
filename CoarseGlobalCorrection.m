function [ output ] = CoarseGlobalCorrection(rhs,tree,domain)

cs = ASCoarseGlobalPreconditioner(tree,domain,rhs);

z = rhs - LaplacianForward(tree,domain,cs);

output = cs + ASPreconditioner(tree,domain,z);
end

