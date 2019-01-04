function [F] = RegCavityFlow(Re,y,leaf)
    
    degs = leaf.degs;
    
    [~,~,~,~,border] = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);    
        
    bump = @(x)exp(1-((x-0.5)/0.5).^2./(1-((x-0.5)/0.5).^2));
    
    Dx = diffmat(degs(1),1,leaf.domain(1,:));
    Dy = diffmat(degs(2),1,leaf.domain(2,:));
    
    Len = prod(degs);
    
    
    u = zeros(degs);
    v = zeros(degs);
    w = zeros(degs);
    
    u(:) = y(1:Len); %x_velocity
    v(:) = y(Len+(1:Len)); %y_velocity
    w(:) = y(2*Len+(1:Len)); %vorticity

    ux = Dx*u; uxx = Dx*ux; uy = u*Dy'; uyy = uy*Dy';
    vx = Dx*v; vxx = Dx*vx; vy = v*Dy'; vyy = vy*Dy';
    wx = Dx*w; wxx = Dx*wx; wy = w*Dy'; wyy = wy*Dy';
   
    P = leaf.points();
      
    f1 = -(uxx+uyy)-wy;
    f1 = f1(:);
    
    f1(border) = u(border) - bump(P(border,1)).*bump(P(border,2)/2);
    
    f2 = -(vxx+vyy)+wx;
    f2 = f2(:);
    f2(border) = v(border);
    
    
    f3 = -1/Re*(wxx+wyy)+u.*wx+v.*wy;
    f3 = f3(:);
    f3(border) = w(border) + uy(border) - vx(border);
    
    F =[f1;f2;f3];
    
end

