function [ F ] = Bratu(u,Approx,kappa,lambda)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    degs = Approx.degs;
    
    %Determine outer border
    [~,~,~,~,border,~] = FindBoundaryIndex2DSides(Approx);
    
    u = reshape(u,degs);
    
    Dx = diffmat(degs(1),1,Approx.domain(1,:));
    Dy = diffmat(degs(2),1,Approx.domain(2,:));

    ux = Dx*u; uxx = Dx*ux; uy = u*Dy'; uyy = uy*Dy';
    
    F = uxx+uyy+kappa*ux+lambda*exp(u)-u;
    
    F(border) = u(border);
    
    F = F(:);
end
