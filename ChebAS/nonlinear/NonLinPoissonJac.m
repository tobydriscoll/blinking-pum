function [ F ] = NonLinPoissonJac(u,Approx)


    F = -Dx*((1+u.^2).*ux)-((1+u.^2).*ux)*Dy';
    

    
    degs = Approx.degs;
    %Determine outer border
    [~,~,~,~,border,sides] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
        
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
    
    J = -Dx*(diag(1+u.^2)*Dx+diag(ux)*(I+diag(2*u)))-Dy*(diag(1+u.^2)*Dy+diag(uy)*(I+diag(2*u)));
    J(sides{4},:) = I(sides{4},:);
    J(sides{1},:) = Dx(sides{1},:);
    J(sides{2},:) = Dx(sides{2},:);
    J(sides{3},:) = Dy(sides{3},:);
end