domain = [0 1;0 1];
deg_in = [5 5];
cheb_struct.domain = domain;
cheb_struct.degs = [32 32];
cheb_struct.split_flag = [true true];
cheb_struct.cdeg_in = deg_in;
cheb_struct.tol = 1e-4;

%Test with 4 patches
Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);
%Tree.split(2);
%Tree.split(1);
%Tree.split(2);

F = PUchebfun(Tree);


tol_n = [1e-4,1e-3];

parms = [20,-1,.5,0];

setInterpMatrices(F);

f = @(u,leaf) CavityFlow(1,u,leaf);
Jac = @(u,leaf) CavityFlowJacobian(1,u,leaf);
[sol, ~, ~, ~] = nsoldPAR_AS(rand(3*length(F),1),f,Jac,F,tol_n,parms);
