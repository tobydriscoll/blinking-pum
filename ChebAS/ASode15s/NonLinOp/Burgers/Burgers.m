
function R = Burgers(NonlinOp,u,nu)

    [U,Ux,Uxx,Uy,Uyy] = NonlinOp.computeDerives(u);

    R = nu*(Uxx+Uyy)-U.*(Ux+Uy);

    R = R(:);
    
end