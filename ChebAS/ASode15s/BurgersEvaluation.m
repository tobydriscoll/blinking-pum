function [ output ] = BurgersEvaluation(Approx,t,y,R)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 
% define f,g
%
 
f = @(x,y,t) 3/4 - 1./(4*(1+exp((R*(4*y-4*x-t)/32))));
g = @(x,y,t) 3/4 + 1./(4*(1+exp((R*(4*y-4*x-t)/32))));
 
    degs = Approx.degs;
    
    %Determine outer border
    [~,out_border,in_border,~] = FindBoundaryIndex2DSides(degs,Approx.domain,Approx.outerbox);
    
    Len = length(Approx);
    
    u = zeros(degs);
    v = zeros(degs);
    
    u(:) = y(1:Len);
    v(:) = y((Len+1):end);
    
    Dx = diffmat(degs(1),1,Approx.domain(1,:));
    Dy = diffmat(degs(2),1,Approx.domain(2,:));
 
    ux = Dx*u; uxx = Dx*ux; uy = u*Dy'; uyy = uy*Dy';
    vx = Dx*v; vxx = Dx*vx; vy = v*Dy'; vyy = vy*Dy';
    
    [X,Y] = ndgrid(Approx.leafGrids{:});
    
    F1 = -(u.*ux+v.*uy) + 1/R*(uxx+uyy);
    F1(out_border) = u(out_border) - f(X(out_border),Y(out_border),t);
    
    
    F2 = -(u.*vx+v.*vy) + 1/R*(vxx+vyy);
    F2(out_border) = v(out_border) - g(X(out_border),Y(out_border),t);
    
    output = [F1(:);F2(:)];
    
end