function [ output ] = CoarseCorrection(rhs,tree,domain,Mat)

cs = ASCoarsePreconditioner(tree,domain,rhs,Mat);

z = rhs - LaplacianForward(tree,domain,cs);

output = cs + ASPreconditioner(tree,domain,z);
end

