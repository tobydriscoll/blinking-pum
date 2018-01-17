function [ output ] = ASPreconditioner(PUApprox,rhs)

LEAVES = PUApprox.collectLeaves();

output = zeros(length(PUApprox),1);

step = zeros(length(LEAVES),1);

for k=2:length(LEAVES)
    step(k) = step(k-1) + length(LEAVES{k-1});
end

for k=1:length(LEAVES)
    
    ind_k = step(k)+(1:length(LEAVES{k}));
    
    rhs_k = rhs(ind_k);
    
    output(ind_k) = LEAVES{k}.linOp\rhs_k;
end