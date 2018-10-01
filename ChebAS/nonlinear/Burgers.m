function [ F ] = Burgers(u,Approx,nu)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    degs = Approx.degs;
    
%    west = @(y) -tanh((-0.1-0.4*y)/(2*nu));
%    east = @(y) -tanh((0.1+0.4*y)/(2*nu));
%    south = @(x) -tanh((x+0.02)/(2*nu));
%    north = @(x) -tanh((x-0.02)/(2*nu));
    
%    west = @(y) 1./(1+exp(y/(2*nu)));
%    east = @(y) 1./(1+exp((1+y)/(2*nu)));
%    south = @(x) 1./(1+exp(x/(2*nu)));
%    north = @(x) 1./(1+exp((1+x)/(2*nu)));
    
    %Determine outer border
    [~,~,~,~,border,~] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
    
    u = reshape(u,degs);
    
    Dx = diffmat(degs(1),1,Approx.domain(1,:));
    Dy = diffmat(degs(2),1,Approx.domain(2,:));

    ux = Dx*u; uxx = Dx*ux; uy = u*Dy'; uyy = uy*Dy';
    
    F = nu*(uxx+uyy)-u.*(ux+uy);
    
    P = Approx.points();
    
%     F(sides{1}) = u(sides{1})-west(P(sides{1},2));
%     F(sides{2}) = u(sides{2})-east(P(sides{2},2));
%     F(sides{3}) = u(sides{3})-south(P(sides{3},1));
%     F(sides{4}) = u(sides{4})-north(P(sides{4},1));
    
    F(border) = u(border) - atan((P(border,1)+P(border,2)-1)*100);
    
    F = F(:);
end
