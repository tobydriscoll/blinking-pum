function [ F ] = NonLinPoisson(u,Approx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    degs = Approx.degs;
    
    
    %Determine outer border
    [~,~,~,~,border,sides] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
    
    north = false(degs); north(:,end) = true;
    rest_border = border & ~ north;
    
    
    u = reshape(u,degs);
    
    Dx = diffmat(degs(1),1,Approx.domain(1,:));
    Dy = diffmat(degs(2),1,Approx.domain(2,:));

    ux = Dx*u; uy = u*Dy';
    
    F = -Dx*((1+u.^2).*ux)-((1+u.^2).*ux)*Dy';
    
    P = Approx.points();
    
    F(sides{4}) = u(sides{4})-1;
    F(sides{1}) = ux(sides{1});
    F(sides{2}) = ux(sides{2});
    F(sides{3}) = uy(sides{3});
    
    F = F(:);
end
