function [ J ] = AllenCahnJacobian(t,u,Approx,ep)
    degs = Approx.degs;
    %Determine outer border
    [~,~,~,~,~,border_s] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
        
    dx = diffmat(degs(1),1,Approx.domain(1,:));
    dy = diffmat(degs(2),1,Approx.domain(2,:));
    
    Ix = eye(degs(1));
    Iy = eye(degs(2));
    
    I = eye(prod(degs));
    
    
    Dx = kron(Iy,dx);
    Dxx = kron(Iy,dx^2);
    
    Dy = kron(dy,Ix);
    Dyy = kron(dy^2,Ix);
    
    J = Dxx+Dyy+(I-diag(3*u.^2))/ep^2;
    
    J(border_s{1},:) = Dx(border_s{1},:);
    J(border_s{2},:) = Dx(border_s{2},:);
    
    J(border_s{3},:) = Dy(border_s{3},:);
    J(border_s{4},:) = Dy(border_s{4},:);
    
    
    
end

