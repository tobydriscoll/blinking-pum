function [fv,J] = SimpNonlinear(u,leaf)

    bound = @(x,y) atan(x.^2+y.^2);

    
    degs = leaf.degs;
    
    [~,~,~,~,border] = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);    
    
    Dx = diffmat(degs(1),1,leaf.domain(1,:));
    Dy = diffmat(degs(2),1,leaf.domain(2,:));

    u  = reshape(u,leaf.degs);
    ux = Dx*u; uxx = Dx*ux; uy = u*Dy'; uyy = uy*Dy';
    
    fv = 0.2*(uxx+uyy)-u.^2;
    
    fv = fv(:);
    points = leaf.points;
    x = points(:,1);
    y = points(:,2);
    
    fv(border) = u(border) - bound(x(border),y(border));
    
end



