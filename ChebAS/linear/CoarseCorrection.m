function [ output ] = CoarseCorrection(PUApprox,rhs,Mat)

cs = ASCoarsePreconditioner(PUApprox,rhs,Mat);

z = rhs - ParSchwarzForward(PUApprox,cs);

output = cs + ASPreconditioner(PUApprox,z);
end

