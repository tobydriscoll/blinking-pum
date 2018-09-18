function [ output ] = ParSchwarzForward(PUApprox,L,U,p,sol)

%PUApprox.sample(sol);

output = zeros(size(sol));

step = zeros(length(PUApprox.leafArray),1);

for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end

for k=1:length(PUApprox.leafArray) 
    
    degs = PUApprox.leafArray{k}.degs;
    
    sol_k = sol(step(k)+(1:prod(degs)));
    
    lap = L{k}*U{k}*sol_k(p{k});
    
    [~,~,in_border,~] = FindBoundaryIndex2DSides(degs,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox);
    
    lap(in_border) = lap(in_border)-PUApprox.leafArray{k}.Binterp*sol;
    
    output(step(k)+(1:length(PUApprox.leafArray{k}))) = lap;
    
end

end


