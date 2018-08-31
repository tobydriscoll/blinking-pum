% INPUT:
%      PUApprox: PUApprox approximation        
%         U,L,p: L U decomposition for jacobian.
%             x: solution
%
% NOTE x is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
function [ output ] = ASPreconditionerMultSols(PUApprox,U,L,p,x)

num_sols = length(x)/length(PUApprox);

output = zeros(length(PUApprox),num_sols);

step = zeros(length(PUApprox.leafArray),1);

x = reshape(x,length(PUApprox),num_sols);

%Figure out starting index for each patch
for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end


%Parallel part
for k=1:length(PUApprox.leafArray)
    
    ind_k = step(k)+(1:length(PUApprox.leafArray{k}));
    
    x_k = x(ind_k,:);
    
    x_k = x_k(:);
    
    output(ind_k,:) = reshape(U{k}\(L{k}\x_k(p{k})),length(PUApprox.leafArray{k}),num_sols);
    
end

output = output(:);
