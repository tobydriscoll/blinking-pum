function [ output ] = ASCoarseGlobalPreconditioner(PUApprox,rhs)

LEAVES = PUApprox.collectLeaves();

%initialize solution to zero
PUApprox.sample(rhs);

x = chebpts(33);

%rhsc = PUApprox.evalfZoneGrid({x x});
rhsc = PUApprox.evalfGrid({x x});

rhsc = rhsc(:);

Dxx = kron(eye(33),diffmat(33,2));
Dyy = kron(diffmat(33,2),eye(33));

Lap = Dxx+Dyy;

border = false(33,33);

border([1 end],:)  = true;
border(:,[1 end])  = true;

border = border(:);

E = eye(33^2);

Lap(border,:) = E(border,:);

solc = Lap\rhsc;

CoarsePatch = ChebPatch([-1 1;-1 1],[-1 1;-1 1],[-1 1;-1 1],[5 5]);

CoarsePatch.sample(solc);

for k=1:length(LEAVES)
    G = CoarsePatch.evalfGrid(LEAVES{k}.leafGrids(),1,0);
    LEAVES{k}.sample(G(:));
end

output = PUApprox.Getvalues();

end