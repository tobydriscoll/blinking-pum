function [F] = SteadyStateBurgers(Re,y,leaf,a,x_0,lambda)
    
    degs = leaf.degs;
    
    u_fun = @(x,y) (1/Re)*(-2*(a(2)+y*a(4)+a(5)*2*lambda*cos(y*lambda).*sinh(lambda*(x-x_0))))./(a(1)+x*a(2)+y*a(3)+x.*y*a(4)+2*a(5)*cos(y*lambda).*cosh(lambda*(x-x_0)));

    v_fun = @(x,y) (1/Re)*(-2*(a(3)+x*a(4)-a(5)*2*lambda*sin(y*lambda).*cosh(lambda*(x-x_0))))./(a(1)+x*a(2)+y*a(3)+x.*y*a(4)+2*a(5)*cos(y*lambda).*cosh(lambda*(x-x_0)));
    
    [~,in_border,out_border,~] = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);    
    
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
    
    f1 = -(uxx+uyy)+Re*(u.*ux+v.*uy);
    f1 = f1(:);
    f1(border) = u(border)-u_fun(P(border,1),P(border,2));
    
    f2 = -(vxx+vyy)+Re*(u.*vx+v.*vy);
    f2 = f2(:);
    f2(border) = v(border)-v_fun(P(border,1),P(border,2));
    
    F =[f1;f2];
    
    
end




