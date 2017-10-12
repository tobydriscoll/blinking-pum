function [ output ] = ASCoarsePreconditioner(PUApprox,domain,rhs,Mat)

LEAVES = PUApprox.collectLeaves({});

%initialize solution to zero
PUApprox.sample(rhs);

PUApprox.Coarsen();

rhsc = PUApprox.Getvalues();

PUApprox.sample(@(x)zeros(length(x),1));

output = [];

solc = Mat\rhsc;

PUApprox.sample(solc);

PUApprox.Refine();

output = PUApprox.Getvalues();

end