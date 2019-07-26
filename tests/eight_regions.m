addpath ..
addpath ../PUChebfun

param.domain = [-1 1;-1 1];
param.degs = [20 20];
param.cdegs = [8 8];
param.split_flag = [true true];
param.tol = 1e-4;
param.odetol = 1e-4;

tspan = [0 5.258];
param.percentClosed = 0.7;
param.pA = 2.14e-6;
param.pS = 6.92e-5;
param.h_e = 2;
param.BoundaryH = 13;
param.initvolume = 34;
param.fluxvolume = 0;

%Test with 16 patches
Tree = ChebPatch(param);
Tree = Tree.split(2);
Tree.children{1} = split(Tree.children{1},1);
Tree.children{2} = split(Tree.children{2},1);
Tree.children{1}.children{1} = split(Tree.children{1}.children{1},2);
Tree.children{1}.children{2} = split(Tree.children{1}.children{2},2);
Tree.children{2}.children{1} = split(Tree.children{2}.children{1},2);
Tree.children{2}.children{2} = split(Tree.children{2}.children{2},2);

Tree.clean();

H = PUchebfun(Tree);

H.sample(@(x,y) zeros(size(x)));
P = H.copy();

setInterpMatrices(H,false);
setInterpMatrices(P,false);

%load initcond_pcl7.mat
[Blinks,M,y0,dy0] = initialize(H,P,param,finalstate);

%%
opt = odeset('mass',M,'reltol',param.odetol,'abstol',param.odetol,...
    'initialstep',1e-5,'initialslope',dy0);
tic
sol = ASode15s(false,Blinks,tspan,y0,{H,P},1,opt);
toc

%%
[h,p,dh,dp] = evaluate(sol,sol.x(end),H,P);
finalstate = struct('H',h,'P',p,'dH',dh,'dP',dp);

