function [ J ] = BurgersJacobian(u,Approx,nu)
    degs = Approx.degs;
    %Determine outer border
    [~,~,~,~,border,~] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
        
    dx = diffmat(degs(1),1,Approx.domain(1,:));
    dy = diffmat(degs(2),1,Approx.domain(2,:));
    
    Ix = eye(degs(1));
    Iy = eye(degs(2));
    
    I = eye(prod(degs));

    Dx = kron(Iy,dx);    
    Dy = kron(dy,Ix);
    
    Dxx = kron(Iy,dx^2);    
    Dyy = kron(dy^2,Ix);
    
    u = reshape(u,degs);
    ux = dx*u; uy = u*dy';
    
    u = u(:); ux = ux(:); uy = uy(:);
    
    J = nu*(Dxx+Dyy)-diag(u)*(Dx+Dy)-diag(ux+uy);
    J(border,:) = I(border,:);
    
end