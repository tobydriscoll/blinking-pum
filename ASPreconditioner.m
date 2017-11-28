function [ output ] = ASPreconditioner(PUApprox,rhs)

LEAVES = PUApprox.collectLeaves({});

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

% PUApprox.sample(output);
% 
% new_output = [];
% for i=1:length(LEAVES)
%     G = PUApprox.evalfGrid(LEAVES{1}.leafGrids(),1,0);
%     new_output = [new_output;G(:)];
% end
% 
% output = new_output;