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
Tree.split(2);

F = PUchebfun(Tree);


tol_n = [1e-4,1e-3];

parms = [20,-1,.5,0];

setInterpMatrices(F);

%[sol, ~, ~, ~] = nsoldPAR_AS(zeros(length(Tree),1),@SimpNonlinear,@SimpNonlinearJac,F,tol_n,parms);


f = @(u,leaf) CavityFlow(1,u,leaf);
Jac = @(u,leaf) CavityFlowJacobian(1,u,leaf);

[sol, ~, ~, ~] = nsoldPAR_AS(zeros(3*length(Tree),1),f,Jac,F,tol_n,parms);

%tol = [1e-10,1e-9];  param = [40, 1000, .5, 0];  % give exact Jacobian
%[u,ithist] = nsold(zeros(length(Tree)*2,1),@(x) Burgers2D(Tree,x,4),tol,param);

%options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',10000,'FunctionTolerance',1e-4);
%[s,~,~,~,~] = fsolve(@(x) SimpNonlinear(Tree,x),zeros(length(Tree),1),options);