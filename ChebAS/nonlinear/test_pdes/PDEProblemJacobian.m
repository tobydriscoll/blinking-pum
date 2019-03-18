function [ J ] = PDEProblemJacobian(u,Approx)
    degs = Approx.degs;
    %Determine outer border
    [~,~,~,~,border,~] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
        
    dx = diffmat(degs(1),1,Approx.domain(1,:));
    dy = diffmat(degs(2),1,Approx.domain(2,:));
    
    Ix = eye(degs(1));
    Iy = eye(degs(2));
    
    I = eye(prod(degs));

    Dxx = kron(Iy,dx^2);    
    Dyy = kron(dy^2,Ix);
    
    J = Dxx+Dyy+diag(3*u.^2);
    
    J(border,:) = I(border,:);
    
end