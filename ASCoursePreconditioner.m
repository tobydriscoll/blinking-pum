function [ output ] = ASCoursePreconditioner(PUApprox,domain,rhs)

LEAVES = PUApprox.collectLeaves({});

%initialize solution to zero
PUApprox.sample(rhs);

PUApprox.Coarsen();

rhsc = PUApprox.Getvalues();

PUApprox.sample(@(x)zeros(length(x),1));

output = [];

step_n = 0;
    
for k=1:length(LEAVES)
    
    dim = LEAVES{k}.degs;
        
    rhs_k = rhsc(step_n+(1:prod(dim)));
    
    pointsl = LEAVES{k}.points();
    
    [out_border,in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
    
    approx = PUApprox.evalf(pointsl(in_border,:),1,0);
    
    bk = [rhs_k(~in_border & ~out_border);approx;rhs_k(out_border)];
    
    step_n = step_n + prod(dim);
    
    LEAVES{k}.sample((LEAVES{k}.ClinOp\bk));
end

step_n = step_n - prod(dim);

for k=length(LEAVES)-1:-1:1
    
    dim = LEAVES{k}.degs;
        
    rhs_k = rhsc(step_n+(1:prod(dim)));
    
    pointsl = LEAVES{k}.points();
    
    [out_border,in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
    
    approx = PUApprox.evalf(pointsl(in_border,:),1,0);
    
    bk = [rhs_k(~in_border & ~out_border);approx;rhs_k(out_border)];
    
    step_n = step_n - prod(dim);
    
    LEAVES{k}.sample((LEAVES{k}.ClinOp\bk));
end

PUApprox.Refine();

output = PUApprox.Getvalues();

end