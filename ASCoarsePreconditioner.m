function [ output ] = ASCoarsePreconditioner(PUApprox,domain,rhs,Mat)

LEAVES = PUApprox.collectLeaves({});

%initialize solution to zero
PUApprox.sample(rhs);

PUApprox.Coarsen();

rhsc = PUApprox.Getvalues();

solc = Mat\rhsc;

PUApprox.sample(solc);

PUApprox.Refine();

output = PUApprox.Getvalues();

end