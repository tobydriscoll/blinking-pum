domain = [-1 1;-1 1];
cheb_struct.domain = domain;
cheb_struct.degs = [15 15];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

pctClosed = 0.78;

BoundaryH = 13;

%Test with 2 patches
 Tree = ChebPatch(cheb_struct);
 Tree = Tree.split(1);
 Tree.split(2);

H = PUchebfun(Tree);
H.sample(@(x,y) zeros(size(x)));

P = H.copy();

setInterpMatrices(H,false);
setInterpMatrices(P,false);

% This sets the blink objects for the patches.
%
% Note in blink, the domain is changed from
% [-1 1]^2 to the local domain.
%
% disc.boundary.all is set to the boundary
% minus the boundary interior to the domain. Maybe
% this is messing things up?
%
% These are the only things that are changed here.
[Blinks,M,y0] = SetBlinks(H,P,pctClosed,BoundaryH);

dt = 0.1;
pred = y0+ParLocalResidual(0,y0,dt,{H,P},Blinks);
%pred = y0;

[J2,L,U,p] = ComputeJacsTime(dt,pred,{H,P},Blinks,dt,M);

E = eye(length(y0));

J = E;

for i=1:length(y0)
    J(:,i) = LinearResidual({H,P},J2,E(:,i));
end

f = @(u) ParLocalResidual(dt,u,dt,{H,P},Blinks)-Masstimes({H,P},M,u-y0);
AJ = jacobi(f,pred);
