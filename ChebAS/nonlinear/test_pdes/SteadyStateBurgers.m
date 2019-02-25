function [F] = SteadyStateBurgers(y,leaf,R)
    
    degs = leaf.degs;
    
    ub = @(x,y) 3/4 - 1./(4*( 1+ exp(R*(-0-4*x+4*y)/32)));
    vb = @(x,y) 3/4 + 1./(4*( 1+ exp(R*(-0-4*x+4*y)/32)));
    
    [~,in_border,out_border,~] = FindBoundaryIndex2DSides(leaf);    
    
    border = in_border | out_border;
    
    Dx = diffmat(degs(1),1,leaf.domain(1,:));
    Dy = diffmat(degs(2),1,leaf.domain(2,:));
    
    Len = prod(degs);
    
    u = zeros(degs);
    v = zeros(degs);
    
    u(:) = y(1:Len); %x_velocity
    v(:) = y(Len+(1:Len)); %y_velocity

    ux = Dx*u; uxx = Dx*ux; uy = u*Dy'; uyy = uy*Dy';
    vx = Dx*v; vxx = Dx*vx; vy = v*Dy'; vyy = vy*Dy';
   
    P = leaf.points();
    
    f1 = 1/R*(uxx+uyy)-(u.*ux+v.*uy)-u;
    f1 = f1(:);
    f1(border) = u(border)-ub(P(border,1),P(border,2));
    
    f2 = 1/R*(vxx+vyy)-(u.*vx+v.*vy)-v;
    f2 = f2(:);
    f2(border) = v(border)-vb(P(border,1),P(border,2));
    
    F =[f1;f2];
    
    
end




