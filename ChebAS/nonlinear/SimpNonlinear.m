function [fv,J] = SimpNonlinear(leaf,u)

    bound = @(x,y) atan(x.^2+y.^2);

    
    degs = leaf.degs;
    
    [~,out_border,in_border,~] = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);    
    
    Dx = kron(eye(degs(2)),diffmat(degs(1),1,leaf.domain(1,:)));
    Dy = kron(diffmat(degs(2),1,leaf.domain(2,:)),eye(degs(1)));
    
    Dxx = kron(eye(degs(2)),diffmat(degs(1),2,leaf.domain(1,:)));
    Dyy = kron(diffmat(degs(2),2,leaf.domain(2,:)),eye(degs(1)));
    
    fv = Dxx*u+Dyy*u-u.^2;
    
    points = leaf.points;
    x = points(:,1);
    y = points(:,2);
    
    fv(out_border) = u(out_border) - bound(x(out_border),y(out_border));
    
    fv(in_border) = u(in_border);
    
    E = eye(prod(degs));
    
    J = Dxx+Dyy-diag(2*u);
    
    J(in_border | out_border,:) = E(in_border | out_border,:);
end



