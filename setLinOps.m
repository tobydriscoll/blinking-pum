function [rhs] = setLinOps(Tree,L,B,force,border)

if Tree.is_leaf
    LEAVES = {Tree};
else
    LEAVES = Tree.collectLeaves({});    
end

step = 0;

rhs = [];

for k=1:length(LEAVES)
    points = LEAVES{k}.points;
    sol = force(points);
    dim = LEAVES{k}.degs;
    
    [out_border_c,out_border,in_border] = FindBoundaryIndex2DSides(dim,LEAVES{k}.domain(),LEAVES{k}.outerbox);
    
    sol(out_border) = border(points(out_border,:));
    sol(in_border) = 0;
    
    rhs = [rhs;sol];
    
    Dx = kron(eye(dim(2)),diffmat(dim(1),1,LEAVES{k}.domain(1,:)));
    Dy = kron(diffmat(dim(2),1,LEAVES{k}.domain(2,:)),eye(dim(1)));
    
    Dxx = kron(eye(dim(2)),diffmat(dim(1),2,LEAVES{k}.domain(1,:)));
    Dyy = kron(diffmat(dim(2),2,LEAVES{k}.domain(2,:)),eye(dim(1)));
    
    E = eye(prod(dim));

    OP = L(E,diag(points(:,1)),diag(points(:,2)),Dx,Dy,Dxx,Dyy);
    
    for i=1:4
        if any(out_border_c{i})
            OP(out_border_c{i},:) = ...
            B{i}(E(out_border_c{i},:),diag(points(out_border_c{i},1)),diag(points(out_border_c{i},2)),Dx(out_border_c{i},:),Dy(out_border_c{i},:),Dxx(out_border_c{i},:),Dyy(out_border_c{i},:));
        end
    end
    
    OP(in_border,:) = E(in_border,:);
    
    LEAVES{k}.linOp = OP;
end

end

