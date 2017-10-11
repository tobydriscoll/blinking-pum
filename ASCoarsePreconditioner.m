function [ output ] = ASCoarsePreconditioner(PUApprox,domain,rhs)

LEAVES = PUApprox.collectLeaves({});

%initialize solution to zero
PUApprox.sample(rhs);

PUApprox.Coarsen();

rhsc = PUApprox.Getvalues();

PUApprox.sample(@(x)zeros(length(x),1));

output = [];



for j=1:4*length(LEAVES)
    
    step_n = 0;
    
    for k=1:length(LEAVES)
        
        dim = LEAVES{k}.degs;
        
        rhs_k = rhsc(step_n+(1:prod(dim)));
        
        pointsl = LEAVES{k}.points();
        
        [out_border,in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
        
        approx = PUApprox.evalfZone(pointsl(in_border,:));
        
        bk = zeros(length(LEAVES{k}),1);
        
        bk(~in_border) = rhs_k(~in_border);
        
        bk(in_border) = approx;
        
        step_n = step_n + prod(dim);
        
        LEAVES{k}.sample((LEAVES{k}.ClinOp\bk));
    end
end

PUApprox.Refine();

output = PUApprox.Getvalues();

end