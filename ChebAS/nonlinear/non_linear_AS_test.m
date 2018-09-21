domain = [0 1;0 1];
deg_in = [5 5];
cheb_struct.domain = domain;
cheb_struct.degs = [50 50];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.cdeg_in = deg_in;
cheb_struct.tol = 1e-4;

tol_n = [1e-5,1e-4];
parms = [20,-1,.5,0];

% %Test with 4 patches
 Tree = ChebPatch(cheb_struct);
 Tree = Tree.split(1);
 Tree.split(2);
 %Tree.split(1);
 %Tree.split(2);
 %Tree.split(1);
 %Tree.split(2);

%Tree = ChebPatch(cheb_struct);
%Tree = Tree.split(2);
%Tree.children{2} = Tree.children{2}.split(2);
%Tree.children{2}.children{2} = Tree.children{2}.children{2}.split(2);
%Tree.children{2}.children{2}.children{2} = Tree.children{2}.children{2}.children{2}.split(2);
%Tree.children{2}.children{2}.children{2}.children{2} = Tree.children{2}.children{2}.children{2}.children{2}.split(2);
%Tree.clean();

leaf_struct.domain = domain;
leaf_struct.degs = [33 33];
leaf_struct.split_flag = [true true];
leaf_struct.cdeg_in = deg_in;
leaf_struct.tol = 1e-4;

Leaf = ChebPatch(leaf_struct);
G = Leaf.leafGrids();

F = PUchebfun(Tree);
F.sample(@(x,y) zeros(size(x)));

setInterpMatrices(F,true);

Re  = 1000;
steep = 0.08;
f = @(u,leaf) CavityFlow(Re,u,leaf,steep);
Jac = @(u,leaf) CavityFlowJacobian(Re,u,leaf);

%f = @(u,leaf) AllenCahn(leaf,0,u,ep);
%Jac = @(u,leaf) AllenCahnJacobian(0,u,leaf,ep);

lambda = 6.808;

f = @(u,leaf) LGB(u,leaf,lambda);
Jac = @(u,leaf) LGBJacobian(u,leaf,lambda);


%f = @ SimpNonlinear;
%Jac = @ SimpNonlinearJac;

% Re = 50;
% a = [2.43e+05 0 0 1.6e5 1];
% lambda = 6;
% x_0 = 1;
%f = @(u,leaf) SteadyStateBurgers(Re,u,leaf,a,x_0,lambda);
%Jac = @(u,leaf) SteadyStateBurgersJacobian(Re,u,leaf);
%init = zeros(length(F)*3,1);


    
init = [];

for i=1:length(F.leafArray)
    
    leaf = F.leafArray{i};
    
    degs = leaf.degs;
    [out_border_s,~,~,~,border] = FindBoundaryIndex2DSides(degs,leaf.domain,leaf.outerbox);
    
    east_west = out_border_s{1} | out_border_s{2};
    south = out_border_s{3};
    north = out_border_s{4};
    
    u = zeros(degs);
    u(north) = 1;
    
    P = leaf.points();
    u = reshape(SideBumpFunc(P(:,2),[0 1],steep),degs);
    
    v = zeros(degs);
    w = zeros(degs);
    
    init = [init;u(:) v(:) w(:)];
    
end

init = init(:,1);

init = zeros(length(F),1);

%F.sample(@(x,y) exp(-((x-0.5)/0.5).^2./(1-((x-0.5)/0.5).^2)).*exp(-((y-0.5)/0.5).^2./(1-((y-0.5)/0.5).^2)));

%init = F.Getvalues();

%

%tic;
%[ sol,normres,normstep,numgm ] = PreconditionedNewtonForward(f,Jac,init,F,[1e-5 1e-5]);
%toc

tic;
[ sol,normres,normstep,numgm ] = PreconditionedNewton(f,Jac,init,F,1e-5,[1e-5 1e-5],Leaf,G);
toc

%tic;
%[ sol,normres,normstep,numgm ] = PreconditionedNewtonTwoLevelG(f,Jac,init,F,[1e-5 1e-5],1e-5,2);
%toc

%tic;
%[ sol,normres,normstep,numgm ] = PreconditionedNewtonTwoLevel3(f,Jac,init,F,Leaf,G,[1e-5 1e-5],1e-5,2);
%toc

%ParPreconditionedTwoLevelG(init,F,f,Jac);

%[sol, ~, ~, ~] = nsoldPAR_AS_two_level(init,f,Jac,F,tol_n,parms);

%tic;
%[sol, it_hist, ierr, x_hist] = nsoldPAR_AS(init,f,Jac,F,[1e-6 1e-5],parms,[1e-10 1e-9]);
%toc

%tic;
%[sol, it_hist, ierr, x_hist] = nsoldPAR_AS_forward(init,f,Jac,F,[1e-10 1e-9],parms);
%toc

%[sol, it_hist, ierr,x_hist] = nsoldAS(init,@(u)f(u,Tree),@(u)Jac(u,Tree),tol_n,parms);

%tol = norm(f(init,Tree))*1e-3;
%options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',5000,'FunctionTolerance',tol,'Display','iter');
%sol = fsolve(@(u)sol_and_jac(@(u)f(u,Tree),@(u)Jac(u,Tree),u),init,options);

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

sol = reshape(sol,length(F),1);
F.sample(sol(:,1));
plot(F);