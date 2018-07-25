function [F] = CavityFlow(Re,y,leaf)
    
    degs = leaf.degs;
    
    [out_border_s,~,~,in_border_s] = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);    
    
    east_west_south = out_border_s{1} | out_border_s{2} | out_border_s{3};
    north = in_border_s{4};
    
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
   
    f1 = -(uxx+uyy)-wy;
    f1 = f1(:);
    f1(east_west_south) = u(east_west_south);
    f1(north) = f1(north) - ones(sum(north),1);
    
    f2 = -(vxx+vyy)+wx;
    f2 = f2(:);
    f2(east_west_south | north) = v(east_west_south | north);
    
    
    f3 = -1/Re*(wxx+wyy)+u.*wx+v.*wy;
    f3 = f3(:);
    f3(east_west_south | north) = f3(east_west_south | north) + uy(east_west_south | north) - vx(east_west_south | north);
    
    F =[f1;f2;f3];
    
    
end



