function [J] = CavityFlowJacobian(Re,y,leaf)
    
    degs = leaf.degs;
    
    Len = prod(degs);
    
    [~,out_border,in_border,~] = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);
    
    border = out_border | in_border;
    
    dx = diffmat(degs(1),1,leaf.domain(1,:));
    dy = diffmat(degs(2),1,leaf.domain(2,:));
    
    Ix = eye(degs(1));
    Iy = eye(degs(2));
    
    w = zeros(degs);
     
    u = y(1:Len); %x_velocity
    v = y(Len+(1:Len)); %y_velocity
    w(:) = y(2*Len+(1:Len)); %vorticity

    wx = dx*w; wy = w*dy';
   
    wx = wx(:); wy = wy(:);
    
    Dx = kron(Iy,dx);
    Dxx = kron(Iy,dx^2);
    
    Dy = kron(dy,Ix);
    Dyy = kron(dy^2,Ix);
    
    I = eye(prod(degs));
    
    Z = zeros(prod(degs));
    
    J1du = -(Dxx+Dyy); J1dv = zeros(Len); J1dw = -Dy;
    J1du(border,:) = I(border,:); J1dv(border,:) = Z(border,:); J1dw(border,:) = Z(border,:); 
    
    J2du = zeros(Len); J2dv = -(Dxx+Dyy); J2dw = Dx;
    J2du(border,:) = Z(border,:); J2dv(border,:) = I(border,:); J2dw(border,:) = Z(border,:); 
    
    J3du = diag(wx);   J3dv = diag(wy);   J3dw = -1/Re*(Dxx+Dyy)+u.*Dx+v.*Dy;
    J3du(border,:) = Dy(border,:); J3dv(border,:) = -Dx(border,:); J3dw(border,:) = I(border,:); 
    
    J = [J1du J1dv J1dw;J2du J2dv J2dw;J3du J3dv J3dw];
    
end

