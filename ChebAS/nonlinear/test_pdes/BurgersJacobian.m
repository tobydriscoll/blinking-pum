function [ J ] = BurgersJacobian(u,Approx,nu)

    %Determine outer border
    [~,~,~,~,border,~] = FindBoundaryIndex2DSides(Approx);
        
    dx = diffmat(Approx.degs(1),1,Approx.domain(1,:));
    dy = diffmat(Approx.degs(2),1,Approx.domain(2,:));
    
    Ix = eye(Approx.degs(1));
    Iy = eye(Approx.degs(2));
    
    I = eye(prod(Approx.degs));

    Dx = kron(Iy,dx);    
    Dy = kron(dy,Ix);
    
    Dxx = kron(Iy,dx^2);    
    Dyy = kron(dy^2,Ix);
    
    u = reshape(u,Approx.degs);
    ux = dx*u; uy = u*dy';
    
    u = u(:); ux = ux(:); uy = uy(:);
    
    J = nu*(Dxx+Dyy)-diag(u)*(Dx+Dy)-diag(ux+uy);
    J(border,:) = I(border,:);
    
end