function [ output ] = CoarseCorrection(rhs,tree,domain)

cs = ASCoarsePreconditioner(tree,domain,rhs);

z = rhs - LaplacianForward(tree,domain,cs);

output = cs + ASPreconditioner(tree,domain,z);
end

