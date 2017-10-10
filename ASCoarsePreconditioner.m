function [ output ] = ASCoarsePreconditioner(PUApprox,domain,rhs)

LEAVES = PUApprox.collectLeaves({});

%initialize solution to zero
PUApprox.sample(rhs);

x = chebpts(33);

rhsc = PUApprox.evalfZoneGrid({x x});

rhsc = rhsc(:);

Dxx = kron(eye(33),diffmat(33,2));
Dyy = kron(diffmat(33,2),eye(33));

Lap = Dxx+Dyy;

border = false(33,33);

border([1 end],:)  = true;
border(:,[1 end])  = true;

border = border(:);

E = eye(33^2);

Lap(border,:) = E(border,:);

solc = Lap\rhsc;

CoarsePatch = ChebPatch([-1 1;-1 1],[-1 1;-1 1],[-1 1;-1 1],[5 5]);

CoarsePatch.sample(solc);

for k=1:length(LEAVES)
    G = CoarsePatch.evalfGrid(LEAVES{k}.leafGrids(),1,0);
    LEAVES{k}.sample(G(:));
end


% PUApprox.Coarsen();
% 
% rhsc = PUApprox.Getvalues();
% 
% PUApprox.sample(@(x)zeros(length(x),1));
% 
% output = [];
% 
% 
% 
% for j=1:1
%     
%     step_n = 0;
%     
%     for k=1:length(LEAVES)
%         
%         dim = LEAVES{k}.degs;
%         
%         rhs_k = rhsc(step_n+(1:prod(dim)));
%         
%         pointsl = LEAVES{k}.points();
%         
%         [out_border,in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
%         
%         approx = PUApprox.evalfZone(pointsl(in_border,:));
%         
%         bk = zeros(length(LEAVES{k}),1);
%         
%         bk(~in_border) = rhs_k(~in_border);
%         
%         bk(in_border) = approx;
%         
%         step_n = step_n + prod(dim);
%         
%         LEAVES{k}.sample((LEAVES{k}.ClinOp\bk));
%     end
%     
% end
% 
% PUApprox.Refine();

output = PUApprox.Getvalues();

end