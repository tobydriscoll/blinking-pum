function [J] = SimpNonlinearJac(u,leaf)

    bound = @(x,y) atan(x.^2+y.^2);
    
    degs = leaf.degs;
    
    Dxx = kron(eye(degs(2)),diffmat(degs(1),2,leaf.domain(1,:)));
    Dyy = kron(diffmat(degs(2),2,leaf.domain(2,:)),eye(degs(1)));
    
    [~,out_border,in_border,~] = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);    
    
    E = eye(prod(degs));
    
    J = 0.2*(Dxx+Dyy)-diag(2*u);
    
    J(in_border | out_border,:) = E(in_border | out_border,:);
end
