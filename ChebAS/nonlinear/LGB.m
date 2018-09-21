function [ F ] = LGB(u,Approx,lambda)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    degs = Approx.degs;
    
    %Determine outer border
    [~,~,~,~,border,~] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
    
    u = reshape(u,degs);
    
    Dx = diffmat(degs(1),1,Approx.domain(1,:));
    Dy = diffmat(degs(2),1,Approx.domain(2,:));

    uxx = Dx^2*u; uyy = u*(Dy')^2;
    
    F = uxx+uyy+lambda*exp(u);
    
    F(border) = u(border);
    
    F = F(:);
end

