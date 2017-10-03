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
    
    Dxx = kron(eye(dim(2)),diffmat(dim(1),2));
    Dyy = kron(diffmat(dim(2),2),eye(dim(1)));
        
    scalex = 2/diff(LEAVES{k}.domain(1,:));
    scaley = 2/diff(LEAVES{k}.domain(2,:));
    
    lap = scalex^2*Dxx+scaley^2*Dyy;
    E1 = eye(prod(dim));
    
    approx = PUApprox.evalf(pointsl(in_border,:),1,0);
    
    Ak = [lap(~(in_border | out_border),:);E1(in_border,:);E1(out_border,:)];
    bk = [rhs_k(~in_border & ~out_border);approx;rhs_k(out_border)];
    
    step_n = step_n + prod(dim);
    
    LEAVES{k}.sample((Ak\bk));
end

PUApprox.Refine();

output = PUApprox.Getvalues();

end