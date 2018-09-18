function [ output ] = ASPreconditioner(PUApprox,L,U,p,rhs)


output = zeros(length(PUApprox),1);

step = zeros(length(PUApprox.leafArray),1);

for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end

for k=1:length(PUApprox.leafArray)
    
    ind_k = step(k)+(1:length(PUApprox.leafArray{k}));
    
    rhs_k = rhs(ind_k);
    
    %output(ind_k) = PUApprox.leafArray{k}.linOp\rhs_k;
    output(ind_k) = U{k}\(L{k}\rhs_k(p{k}));
end