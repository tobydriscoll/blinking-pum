domain = [-1 1;-1 1];
cheb_struct.domain = domain;
cheb_struct.degs = [33 33];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

%Test with 16 patches
%  Tree = ChebPatch(cheb_struct);
%  Tree = Tree.split(1);
%  Tree.split(2);
%  Tree.split(1);
%  Tree.split(2);
%  Tree.split(1);
%  Tree.split(2);
%  
% F = PUchebfun(Tree);
% F.sample(@(x,y) zeros(size(x)));

%F = PUchebfun(@(x,y)atan((x+y)*5),[-1 1;-1 1],'Degree',[33 33],'tol',1e-7);
%F = PUchebfun(@(x,y)SideBumpFunc(y,[0 1],0.1),[0 1;0 1],'Degree',[33 33],'CoarseDegree',[9 9],'tol',1e-7);
F = PUchebfun(@(x,y)exp(-x.^20./(1-x.^20)).*exp(-y.^20./(1-y.^20)),[-1 1;-1 1],'Degree',[33 33],'tol',5e-5);

F.reset;

setInterpMatrices(F,true);

% Re  = 1000;
% steep = 0.1;
% f = @(u,leaf) CavityFlow(Re,u,leaf,steep);
% Jac = @(u,leaf) CavityFlowJacobian(Re,u,leaf);
% 
% F.Setvalues(@(x,y)SideBumpFunc(y,[0 1],steep)); %set values
% F.sample(@(x,y)SideBumpFunc(y,[0 1],steep)); %set coeffs
% u = F.Getvalues();
% Fy = diff(F,2,1);
% w = -(Fy.Getvalues());
% v = zeros(length(F),1);
% init = [u;v;w];


% bound_f = @(x,y) atan((cos(pi*3/16)*x+sin(pi*3/16)*y)*1);
% nu = 1/2;
% f = @(u,leaf) Burgers(u,leaf,nu,bound_f);
% Jac = @(u,leaf) BurgersJacobian(u,leaf,nu);
% F.Setvalues(bound_f);
% init = F.Getvalues();

bound_f = @(x,y) atan((cos(pi*3/16)*x+sin(pi*3/16)*y)*1);
nu = 1/100;
f = @(u,leaf) Burgers(u,leaf,nu,bound_f);
Jac = @(u,leaf) BurgersJacobian(u,leaf,nu);
F.sample(bound_f);
init = F.Getvalues();


%  tic;
%  [ sol,normres1,normstep1,numgm1 ] = NSKsolver(f,Jac,init,F,[1e-10 1e-10]);
%  toc

% %   tic;
%  [ sol,normres2,normstep2,numgm2,normresf2 ] = SNKsolver(f,Jac,init,F,[1e-10 1e-10]);
% toc

tic;
[ sol,normres3,normstep3,numgm3 ] = SNK2levelsolver(f,Jac,init,F,[1e-10 1e-10],1e-3,2);
toc

solr = reshape(sol,length(F),1);
F.sample(solr(:,1));
plot(F);