domain = [0 1;0 1];
cheb_struct.domain = domain;
cheb_struct.degs = [33 33];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

tol_n = [1e-10 1e-10];

Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);
Tree.split(2);
%Tree.split(1);
%Tree.split(2);

F = PUchebfun(Tree);
F.sample(@(x,y) zeros(size(x)));

setInterpMatrices(F,true);

bound_f = @(x,y) atan((cos(pi*3/16)*x+sin(pi*3/16)*y)*1);

nu = 1/1500;

f = @(NonLinOp,u)Burgers(NonLinOp,u,nu);
Jac = @(NonLinOp,u)BurgersJacobian(NonLinOp,u,nu);

NonLinOps = SetUpNonLinOps(F,f,Jac,bound_f);

F.pack();

F.sample(bound_f);
init = F.Getvalues();

tic;
 [ sol,normres2,normstep2,numgm2,tolg2 ] = NSKsolver(init,F,NonLinOps,tol_n);
toc

sol = F.Getunpackedvalues(sol);

F.unpack();

F.sample(sol);
plot(F);