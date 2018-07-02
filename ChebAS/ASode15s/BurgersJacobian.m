function [ J ] = BurgersJacobian(t,y,Approx,R)
    degs = Approx.degs;
    %Determine outer border
    [~,out_border,in_border,~] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
    
    border = out_border | in_border;
    
    Len = length(Approx);
    
    u = y(1:Len);
    v = y((Len+1):end);
    
    dx = diffmat(degs(1),1,Approx.domain(1,:));
    dy = diffmat(degs(2),1,Approx.domain(2,:));
    
    Ix = eye(degs(1));
    Iy = eye(degs(2));
    
    I = eye(prod(degs));
    
    Z = zeros(prod(degs));
    
    Dx = kron(Iy,dx);
    Dxx = kron(Iy,dx^2);
    
    Dy = kron(dy,Ix);
    Dyy = kron(dy^2,Ix);
    
    J1du = -(diag(u)*Dx+diag(Dx*u)+diag(v)*Dy)+1/R*(Dxx+Dyy);
    J1du(border,:) = I(border,:);
    
    J1dv = -diag(Dy*u);
    J1dv(border,:) = Z(border,:);
    
    J2du = -diag(Dx*v);
    J2du(border,:) = Z(border,:);
    
    J2dv = -(diag(u)*Dx+diag(v)*Dy+diag(Dy*v))+1/R*(Dxx+Dyy);
    J2dv(border,:) = I(border,:);
    
    J = [J1du J1dv;J2du J2dv];
    
end