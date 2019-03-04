domain = [-1 1;-1 1];
cheb_struct.domain = domain;
cheb_struct.degs = [33 33];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

pctClosed = 0.75;
pA = 2.14e-2;
pS = 6.92e-4;
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
[Blinks,M,y0] = setBlinks(H,P,pctClosed,BoundaryH,pA,pS);

%This move [H_1 H_2 P_1 P_2] to
%{[H_1 P_1],[H_2 P_2]}
[ sol_loc ] = unpackPUvecs(y0,{H P});

%Test on first block
t = 0.0377;
residual = @(y) Blinks{1}.timederiv(t,y);
Jac = @(y) Blinks{1}.jac(t,y);

AJ = jacobi(residual,sol_loc{1});
J = Jac(sol_loc{1});

y02 = packPUvecs(sol_loc,{H P});

