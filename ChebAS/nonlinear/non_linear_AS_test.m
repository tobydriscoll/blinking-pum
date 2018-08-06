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
bound = @(x,y) atan(x.^2+y.^2);

F = PUchebfun(Tree);


tol_n = [1e-5,1e-4];

parms = [20,-1,.5,0];

F.sample(@(x,y) zeros(size(x)));

setInterpMatrices(F,true);

%f = @(u,leaf) CavityFlow(1,u,leaf);
%Jac = @(u,leaf) CavityFlowJacobian(1,u,leaf);

f = @ SimpNonlinear;
Jac = @ SimpNonlinearJac;

%[sol, ~, ~, ~] = nsoldPAR_AS(rand(length(F),1),f,Jac,F,tol_n,parms);



%init = rand(length(F),1);

F.sample(bound);
init = F.Getvalues();

[sol, ~, ~, ~] = nsoldPAR_AS_two_level(zeros(length(F),1),f,Jac,F,tol_n,parms);


F.sample(sol);