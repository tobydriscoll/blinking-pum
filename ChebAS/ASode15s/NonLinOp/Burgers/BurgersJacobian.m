function [ J ] = BurgersJacobian(NonlinOp,u,nu)
%BURGERSJACOBIAN Summary of this function goes here
%   Detailed explanation goes here

[U,Ux,~,Uy,~] = NonlinOp.computeDerives(u);

D = NonlinOp.GetDisc();

J = nu*(D.Dx2+D.Dy2) - diag(U(:))*(D.Dx+D.Dy)-diag(Ux(:)+Uy(:));

end

