function [J] = SteadyStateBurgersJacobian(Re,y,leaf)

degs = leaf.degs;
    
    Len = prod(degs);
    
    [~,out_border,in_border,~] = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);
    
    border = out_border | in_border;
    
    dx = diffmat(degs(1),1,leaf.domain(1,:));
    dy = diffmat(degs(2),1,leaf.domain(2,:));
    
    Ix = eye(degs(1));
    Iy = eye(degs(2));
    
    u = zeros(degs);
    v = zeros(degs);
    
    u(:) = y(1:Len); %x_velocity
    v(:) = y(Len+(1:Len)); %y_velocity

    ux = dx*u; uy = u*dy';
    vx = dx*v; vy = v*dy';
    
    u = u(:); v = v(:);
    ux = ux(:); uy = uy(:);
    vx = vx(:); vy = vy(:);
    
    
    Dx = kron(Iy,dx);
    Dxx = kron(Iy,dx^2);
    
    Dy = kron(dy,Ix);
    Dyy = kron(dy^2,Ix);
    
    I = eye(prod(degs));
    
    Z = zeros(prod(degs));
    
    J1du = -(Dxx+Dyy)+Re*(diag(u)*Dx+diag(ux)+diag(v)*Dy); J1dv = Re*diag(uy);
    J1du(border,:) = I(border,:); J1dv(border,:) = Z(border,:);
    
    J2du = Re*diag(vx); J2dv = -(Dxx+Dyy)+Re*(diag(u)*Dx+diag(v)*Dy+diag(vy));
    J2du(border,:) = Z(border,:); J2dv(border,:) = I(border,:); 
    
    
    J = [J1du J1dv;J2du J2dv];
end