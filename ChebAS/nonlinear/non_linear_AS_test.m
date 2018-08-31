domain = [-1 1;-0.5 0.5];
deg_in = [5 5];
cheb_struct.domain = domain;
cheb_struct.degs = [33 33];
cheb_struct.split_flag = [true true];
cheb_struct.cdeg_in = deg_in;
cheb_struct.tol = 1e-4;

%Test with 4 patches
Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);
Tree.split(2);
%Tree.split(1);
%Tree.split(2);
%Tree.split(1);
%Tree.split(2);

% Re = 50;
% a = [2.43e+05 0 0 1.6e5 1];
% lambda = 6;
% x_0 = 1;

Re = 10;
a = [20 10 -7 3 5];
lambda = 1;
x_0 = 1;

bound = @(x,y) atan(x.^2+y.^2);

F = PUchebfun(Tree);


tol_n = [1e-10,1e-9];

parms = [20,-1,.5,0];

F.sample(@(x,y) zeros(size(x)));

setInterpMatrices(F,true);

%f = @(u,leaf) CavityFlow(1,u,leaf);
%Jac = @(u,leaf) CavityFlowJacobian(1,u,leaf);

%f = @ SimpNonlinear;
%Jac = @ SimpNonlinearJac;

f = @(u,leaf) SteadyStateBurgers(Re,u,leaf,a,x_0,lambda);
Jac = @(u,leaf) SteadyStateBurgersJacobian(Re,u,leaf);

init = zeros(length(F)*2,1);

%F.sample(bound);
%init = F.Getvalues();


tic;
[ sol,normres,normstep,numgm ] = PreconditionedNewtonForward(f,Jac,init,F,1e-10);
toc

% tic;
% [ sol,normres,normstep,numgm ] = PreconditionedNewton(f,Jac,init,F,1e-10);
% toc

%[ sol,normres,normstep,numgm ] = PreconditionedNewtonTwoLevel(f,Jac,init,F,1e-5,2);

%ParPreconditionedTwoLevelG(init,F,f,Jac);

%[sol, ~, ~, ~] = nsoldPAR_AS_two_level(init,f,Jac,F,tol_n,parms);

% for i=1:1000000
%     new_sol = ParNonLinForward(init,F,f,Jac);
%     res(i) = norm(new_sol-init);
%     init = new_sol;
%     
%     if res(i)<1e-10
%         break;
%     end
%     
%     res(i)
% end

sol = reshape(sol,length(F),2);
F.sample(sol(:,1));
plot(F);