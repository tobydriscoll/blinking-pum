function [ output ] = ASPreconditioner(PUApprox,domain,rhs)

LEAVES = PUApprox.collectLeaves({});

output = [];

step_n = 0;

for k=1:length(LEAVES)
    
    dim = LEAVES{k}.degs;
        
    rhs_k = rhs(step_n+(1:prod(dim)));
    
    step_n = step_n + prod(dim);
    
    output = [output;LEAVES{k}.linOp\rhs_k];
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