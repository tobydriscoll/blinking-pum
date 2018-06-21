function [ output ] = ASCoarsePreconditioner(PUApprox,rhs,Mat)

%initialize solution to zero
PUApprox.Setvalues(rhs);

PUApprox.Coarsen();

rhsc = PUApprox.Getvalues();

solc = Mat\rhsc;

PUApprox.Setvalues(solc);

PUApprox.Refine();

output = PUApprox.Getvalues();

end