domain = [0 1;0 1];
cheb_struct.domain = domain;
cheb_struct.degs = [33 33];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.cdeg_in = deg_in;
cheb_struct.tol = 1e-4;

tol_n = [1e-5,1e-4];
parms = [1000,-1,.5,0];

% %Test with 16 patches
 Tree = ChebPatch(cheb_struct);
 Tree = Tree.split(1);
 Tree.split(2);
 Tree.split(1);
 Tree.split(2);
 
% Tree = ChebPatch(cheb_struct);
% Tree = Tree.split(2);
% Tree.children{2} = Tree.children{2}.split(2);
% Tree.children{2}.children{2} = Tree.children{2}.children{2}.split(2);
% Tree.children{2}.children{2}.children{2} = Tree.children{2}.children{2}.children{2}.split(2);
% Tree.children{2}.children{2}.children{2}.children{2} = Tree.children{2}.children{2}.children{2}.children{2}.split(2);
Tree.clean();

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

Re  = 10;
steep = 0.1;
f = @(u,leaf) CavityFlow(Re,u,leaf,steep);
Jac = @(u,leaf) CavityFlowJacobian(Re,u,leaf);

F.Setvalues(@(x,y)SideBumpFunc(y,[0 1],steep)); %set values
F.sample(@(x,y)SideBumpFunc(y,[0 1],steep)); %set coeffs
u = F.Getvalues();
Fy = diff(F,2,1);
w = -(Fy.Getvalues());
v = zeros(length(F),1);

init = [u;v;w];

% tic;
% [ sol,normres1,normstep1,numgm1 ] = NSKsolver(f,Jac,init,F,[1e-7 1e-7]);
% toc

%  tic;
% [ sol,normres2,normstep2,numgm2,normresf2 ] = SNKsolver(f,Jac,init,F,[1e-7 1e-7]);
% toc
% 
tic;
[ sol,normres3,normstep3,numgm3 ] = SNK2levelsolver(f,Jac,init,F,[1e-10 1e-10],1e-3,2);
toc

solr = reshape(sol,length(F),3);
F.sample(solr(:,1));
plot(F);