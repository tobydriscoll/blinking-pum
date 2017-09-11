function [ output ] = RASPreconditioner(PUApprox,domain,rhs)

LEAVES = PUApprox.collectLeaves({});

output = [];

step_n = 0;

for k=1:length(LEAVES)
    
    dim = LEAVES{k}.degs;
        
    rhs_k = rhs(step_n+(1:prod(dim)));
        
    [out_border,in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
    
    Dxx = kron(eye(dim(2)),diffmat(dim(1),2));
    Dyy = kron(diffmat(dim(2),2),eye(dim(1)));
        
    scalex = 2/diff(LEAVES{k}.domain(1,:));
    scaley = 2/diff(LEAVES{k}.domain(2,:));
    
    lap = scalex^2*Dxx+scaley^2*Dyy;
    E1 = eye(prod(dim));
    A1 = [lap(~(in_border | out_border),:);E1(in_border,:);E1(out_border,:)];
    b1 = [rhs_k(~in_border & ~out_border);zeros(sum(in_border),1);rhs_k(out_border)];
    
    step_n = step_n + prod(dim);
    
    output = [output;(A1\b1)];
end

end