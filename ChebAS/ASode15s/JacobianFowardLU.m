% INPUT:
%      PUApprox: PUApprox approximation        
%         J: cell array of local Jacobians
%         x: solution
function [ output ] = JacobianFowardLU(PUApprox,L,U,p,x,num_sols)


output = zeros(length(PUApprox)*num_sols,1);

step = zeros(length(PUApprox.leafArray),1);

%Figure out starting index for each patch
for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + num_sols*length(PUApprox.leafArray{k-1});
end

%non parallel part
for k=1:length(PUApprox.leafArray) 
    degs = PUApprox.leafArray{k}.degs;
    [~,~,in_border,~] = FindBoundaryIndex2DSides(degs,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox);
    
    b = [];
    
    sol_len = length(PUApprox);
    
    for i=1:num_sols
        b_loc = zeros(length(PUApprox.leafArray{k}),1);
        
        b_loc(in_border) = PUApprox.leafArray{k}.Binterp*x((i-1)*sol_len+(1:sol_len));
        b = [b;b_loc];
    end
    
    z{k} = b;
end

%Parallel part
for k=1:length(PUApprox.leafArray)
    
    ind_k = step(k)+(1:num_sols*length(PUApprox.leafArray{k}));
    
    x_k = x(ind_k);
    
    output(ind_k) = U{k}\(L{k}\z{k}(p{k}))-x_k;
    
end

