function [ output ] = RASPreconditioner(PUApprox,domain,force,border,sol)

PUApprox.sample(sol);

LEAVES = PUApprox.collectLeaves({});

output = [];

for k=1:length(LEAVES)
    
    points = LEAVES{k}.points;
    
    dim = LEAVES{k}.degs;
        
    [out_border,in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
    
    Dxx = kron(eye(dim(2)),diffmat(dim(1),2));
    Dyy = kron(diffmat(dim(2),2),eye(dim(1)));
    
    bval = PUApprox.evalf(points(in_border,:),1,0);
    outbval = border(points(out_border,:));
    
    scalex = 2/diff(LEAVES{k}.domain(1,:));
    scaley = 2/diff(LEAVES{k}.domain(2,:));
    
    lap = scalex^2*Dxx+scaley^2*Dyy;
    E1 = eye(prod(dim));
    A1 = [lap(~(in_border | out_border),:);E1(in_border,:);E1(out_border,:)];
    b1 = [force(points(~in_border & ~out_border,:));bval;outbval];
    
    output = [output;(A1\b1)];
end

end

