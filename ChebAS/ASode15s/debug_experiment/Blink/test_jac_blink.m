domain = [-1 1;-1 1];
cheb_struct.domain = domain;
cheb_struct.degs = [10 10];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

pctClosed = 0.75;
pA = 0;
pS = 1e-3;
BoundaryH = 13;

%Test with 2 patches
 Tree = ChebPatch(cheb_struct);
 Tree = Tree.split(1);
 %Tree.split(2);

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
[Blinks,M,y0,result] = setBlinks(H,P,pctClosed,BoundaryH,pA,pS,2);

%This move [H_1 H_2 P_1 P_2] to
%{[H_1 P_1],[H_2 P_2]}
[ sol_loc ] = unpackPUvecs(y0,{H P});

%Test on first block
t = 0;
residual = @(y) Blinks{1}.timederiv(t,y);
Jac = @(y) Blinks{1}.jac(t,y);

E = eye(length(y0));

[Jacs,l,u,p] = ComputeJacsTime(0.1,y0,{H,P},Blinks,0.1,M);


JJ = E;

for i=1:length(y0)
   JJ(:,i) = LinearResidual({H,P},Jacs,E(:,i));
end

J = ASJacTime({H,P},Blinks,M,0.1,0.1,y0);
AJ = jacobi(residual,sol_loc{1});
% J = Jac(sol_loc{1});
% 
% ind = false(length(sol_loc{1}),1);
% 
% H_ind = ind;
% 
% H_ind(1:length(H.leafArray{1}))=true;
% 
% P_ind = ~H_ind;

% 
% y02 = packPUvecs(sol_loc,{H P});

