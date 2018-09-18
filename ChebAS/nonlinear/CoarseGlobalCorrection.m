function [ output ] = CoarseGlobalCorrection(rhs,J,L,U,p,Jac_hat,PUApprox,Leaf,G)

cs = LinearGlobalCoarseCorrect( PUApprox, Leaf, G, rhs,Jac_hat);

z = rhs - ParLinearResidual(cs,PUApprox,J);

output = cs + ASPreconditionerMultSols(PUApprox,U,L,p,z);
end