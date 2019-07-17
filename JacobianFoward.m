% INPUT:
%      PUApprox: PUApprox approximation        
%         J: cell array of local Jacobians
%         x: solution
function [ output ] = JacobianFoward(PUApprox,J,x)


num_sols = length(x)/length(PUApprox);

output = zeros(length(PUApprox),num_sols);

step = zeros(length(PUApprox.leafArray),1);

x = reshape(x,length(PUApprox),num_sols);

%Figure out starting index for each patch
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

%Parallel part
for k=1:length(PUApprox.leafArray)
    
    ind_k = step(k)+(1:length(PUApprox.leafArray{k}));
    
    output(ind_k) = J{k}\z{k}-x(ind_k);
    
end

