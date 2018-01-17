function [ output ] = ASCoarsePreconditioner(PUApprox,rhs,Mat)

%initialize solution to zero
PUApprox.sample(rhs);

PUApprox.Coarsen();

rhsc = PUApprox.Getvalues();

solc = Mat\rhsc;

PUApprox.sample(solc);

PUApprox.Refine();

output = PUApprox.Getvalues();

end