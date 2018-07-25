function [fv,J] = SimpNonlinear(u,leaf)

    bound = @(x,y) atan(x.^2+y.^2);

    
    degs = leaf.degs;
    
    [~,out_border,in_border,~] = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);    
    
    Dx = diffmat(degs(1),1,leaf.domain(1,:));
    Dy = diffmat(degs(2),1,leaf.domain(2,:));

    u  = reshape(u,leaf.degs);
    ux = Dx*u; uxx = Dx*ux; uy = u*Dy'; uyy = uy*Dy';
    
    fv = uxx+uyy-u.^2;
    
    fv = fv(:);
    points = leaf.points;
    x = points(:,1);
    y = points(:,2);
    
    fv(out_border) = u(out_border) - bound(x(out_border),y(out_border));
    
    fv(in_border) = u(in_border);
    
end



