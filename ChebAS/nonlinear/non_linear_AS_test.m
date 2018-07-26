domain = [0 1;0 1];
deg_in = [5 5];
cheb_struct.domain = domain;
cheb_struct.degs = [5 5];
cheb_struct.split_flag = [true true];
cheb_struct.cdeg_in = deg_in;
cheb_struct.tol = 1e-4;

%Test with 4 patches
Tree = ChebPatch(cheb_struct);
%Tree = Tree.split(1);
%Tree.split(2);

F = PUchebfun(Tree);


tol_n = [1e-4,1e-3];

parms = [20,-1,.5,0];

setInterpMatrices(F);

%[sol, ~, ~, ~] = nsoldPAR_AS(zeros(length(Tree),1),@SimpNonlinear,@SimpNonlinearJac,F,tol_n,parms);


%f = @(u,leaf) CavityFlow(1,u,leaf);
%Jac = @(u,leaf) CavityFlowJacobian(1,u,leaf);
%[sol, ~, ~, ~] = nsoldPAR_AS(zeros(3*length(F),1),f,Jac,F,tol_n,parms);

tol = [1e-10,1e-9];  param = [40, 1000, .5, 0];  % give exact Jacobian

f = @(u) CavityFlow(1,u,Tree);
Jac = @(u) CavityFlowJacobian(1,u,Tree);
%sol = nsoldAS(zeros(length(F)*3,1),f,Jac,tol,parms);

options = optimoptions('fsolve','SpecifyObjectiveGradient',true);
sol = fsolve(@(u)sol_and_jac(f,Jac,u),rand(length(F)*3,1),options);

%options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',10000,'FunctionTolerance',1e-4);
%[s,~,~,~,~] = fsolve(@(x) SimpNonlinear(Tree,x),zeros(length(Tree),1),options);