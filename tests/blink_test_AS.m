addpath ..
addpath ../PUChebfun

param.domain = [-1 1;-1 1];
param.degs = [25 25];
param.cdegs = [9 9];
param.split_flag = [true true];
param.tol = 1e-4;
param.odetol = 1e-4;

tspan = [0 0.02];
param.percentClosed = 0.7;
model.pA = 2.14e-6;
model.pS = 6.92e-5;
param.h_e = 2;
param.BoundaryH = 13;
param.initvolume = 34;
param.fluxvolume = 0;

%Test with 4 patches
Tree = ChebPatch(param);
overlap = 0.2;

%This replaces Tree with one leaf overlap away from a boundary
%| |      |
%| |      |
%| |      |
%| |      |
Tree = Tree.split(2,false,overlap);

%This does the same but with the larger leaf
%| |     | |
%| |     | |
%| |     | |
%| |     | |
Tree.children{2} = Tree.children{2}.split(2,false,overlap);


%This does the same but with the larger leaf
%| |     | |
%| |     | |
%| |_ _ _| |
%| |     | |
Tree.children{2}.children{1} = Tree.children{2}.children{1}.split(1,false,overlap);


%This does the same but with the larger left over leaf
%| |_ _ _| |
%| |     | |
%| |_ _ _| |
%| |     | |
Tree.children{2}.children{1}.children{2} = Tree.children{2}.children{1}.children{2}.split(1,false,overlap);

Tree.clean();

H = PUchebfun(Tree);

%H = PUchebfun(@(x,y)exp(-x.^20./(1-x.^20)).*exp(-y.^20./(1-y.^20)),[-1 1;-1 1],'Degree',[20 20],'CoarseDegree',[9 9],'tol',1e-3);
%H.reset();

H.sample(@(x,y) zeros(size(x)));
P = H.copy();

setInterpMatrices(H,false);
setInterpMatrices(P,false);

load initcond_pcl7.mat
[Blinks,M,y0,dy0] = initialize(H,P,param);

opt = odeset('mass',M,'reltol',param.odetol,'abstol',param.odetol,...
    'initialstep',1e-10,'initialslope',dy0);

[t,U] = ASode15s(true,Blinks,tspan,y0,{H,P},1,opt);



