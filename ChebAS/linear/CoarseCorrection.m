function [ output ] = CoarseCorrection(PUApprox,rhs,M,L,U,p,Mat)

cs = ASCoarsePreconditioner(PUApprox,rhs,Mat);

z = rhs - ParSchwarzForward(PUApprox,M,cs);

output = cs + ASPreconditioner(PUApprox,L,U,p,z);
end

