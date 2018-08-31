% INPUT:
%      PUApprox: PUApprox approximation        
%         J: cell array of local Jacobians
%         x: solution
%
% NOTE x is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
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
    
    z{k} = zeros(prod(degs),num_sols);
    
    
    if ~PUApprox.iscoarse
        z{k}(in_border,:) = PUApprox.leafArray{k}.Binterp*x;
    else
        z{k}(in_border,:) = PUApprox.leafArray{k}.CBinterp*x;
    end
    
end

%Parallel part
for k=1:length(PUApprox.leafArray)
    
    ind_k = step(k)+(1:length(PUApprox.leafArray{k}));
    
    x_k = x(ind_k,:);
    
   
    output(ind_k,:) =  reshape(J{k}*x_k(:),length(PUApprox.leafArray{k}),num_sols) - z{k};
    
end

output = output(:);


