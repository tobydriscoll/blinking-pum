function [ J ] = BratuJacobian(u,Approx,kappa,lambda)
    degs = Approx.degs;
    %Determine outer border
    [~,~,~,~,border,~] = FindBoundaryIndex2DSides(Approx);
        
    dx = diffmat(degs(1),1,Approx.domain(1,:));
    dy = diffmat(degs(2),1,Approx.domain(2,:));
    
    Ix = eye(degs(1));
    Iy = eye(degs(2));
    
    I = eye(prod(degs));

    Dx = kron(Iy,dx);    
    
    Dxx = kron(Iy,dx^2);    
    Dyy = kron(dy^2,Ix);
    
    J = Dxx+Dyy+kappa*Dx+lambda*diag(exp(u))-eye(length(u));
    
    J(border,:) = I(border,:);
    
end