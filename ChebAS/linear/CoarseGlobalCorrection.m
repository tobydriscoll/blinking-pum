function [ output ] = CoarseGlobalCorrection(rhs,L,U,p,PUApprox)

cs = LinearGlobalCoarseCorrect( PUApprox, Leaf, G, rhs,Jac_hat);

z = rhs - ParSchwarzForward(PUApprox,L,U,p,cs);

output = cs + ASPreconditioner(PUApprox,L,U,p,z);
end

