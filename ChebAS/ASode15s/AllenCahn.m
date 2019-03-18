function [ output ] = AllenCahn(Approx,t,u,ep)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% define f,g
%
    degs = Approx.degs;
    
    %Determine outer border
    [~,~,~,~,border,border_s] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
    
    u = reshape(u,degs);
    
    Dx = diffmat(degs(1),1,Approx.domain(1,:));
    Dy = diffmat(degs(2),1,Approx.domain(2,:));

    ux = Dx*u; uxx = Dx*ux; uy = u*Dy'; uyy = uy*Dy';
    
    output = uxx+uyy+(u-u.^3)/ep^2;
    
    output(border_s{1}) = ux(border_s{1});
    output(border_s{2}) = ux(border_s{2});
    output(border_s{3}) = uy(border_s{3});
    output(border_s{4}) = uy(border_s{4});
    
    output = output(:);
    
end

