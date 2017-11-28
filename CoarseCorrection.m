function [ output ] = CoarseCorrection(rhs,tree,Mat)

cs = ASCoarsePreconditioner(tree,rhs,Mat);

z = rhs - LaplacianForward(tree,cs);

output = cs + ASPreconditioner(tree,z);
end

