function [ F ] = Burgers(u,Approx,nu,bound_f)


    [~,~,~,~,border,~] = FindBoundaryIndex2DSides(Approx);
    
    u = reshape(u,Approx.degs);
    
    Dx = diffmat(Approx.degs(1),1,Approx.domain(1,:));
    Dy = diffmat(Approx.degs(2),1,Approx.domain(2,:));

    ux = Dx*u; uxx = Dx*ux; uy = u*Dy'; uyy = uy*Dy';
    
    F = nu*(uxx+uyy)-u.*(ux+uy);
    
    P = Approx.points();
    
%     F(sides{1}) = u(sides{1})-west(P(sides{1},2));
%     F(sides{2}) = u(sides{2})-east(P(sides{2},2));
%     F(sides{3}) = u(sides{3})-south(P(sides{3},1));
%     F(sides{4}) = u(sides{4})-north(P(sides{4},1));
    
    F(border) = u(border) - bound_f(P(border,1),P(border,2));
    
    F = F(:);
end
