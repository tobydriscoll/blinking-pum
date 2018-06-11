function [ output ] = JacobianFoward(PUApprox,J,x)


output = zeros(length(PUApprox),1);

step = zeros(length(PUApprox.leafArray),1);

for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end

%non parallel part
for k=1:length(PUApprox.leafArray) 
    degs = PUApprox.leafArray{k}.degs;
    [~,~,in_border,~] = FindBoundaryIndex2DSides(degs,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox);
    z{k} = zeros(length(PUApprox.leafArray{k}),1);
    z{k}(in_border) = PUApprox.leafArray{k}.Binterp*x;
end

for k=1:length(PUApprox.leafArray)
    
    ind_k = step(k)+(1:length(PUApprox.leafArray{k}));
    
    output(ind_k) = J{k}\z{k}-x(ind_k);
    
end

