function [ F ] = PDEProblem(u,Approx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    degs = Approx.degs;
    
    %Determine outer border
    [~,~,~,~,border,~] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
    
    u = reshape(u,degs);
    
    Dx = diffmat(degs(1),1,Approx.domain(1,:));
    Dy = diffmat(degs(2),1,Approx.domain(2,:));

    ux = Dx*u; uxx = Dx*ux; uy = u*Dy'; uyy = uy*Dy';
    
    F = uxx+uyy+u.^3;
    
    F(border) = u(border);
    
    F = F(:);
end

