addpath ..
addpath ../PUChebfun

param.domain = [-1 1;-1 1];
param.degs = [15 15];
param.cdegs = [8 8];
param.split_flag = [true true];
param.tol = 1e-4;
param.odetol = 1e-4;

prob = [];
prob.percentclosed = 0;
prob.A = 2e-6;
prob.S = 6.92e-5;
prob.h_slideover = 2;
prob.h_boundary = 13;
prob.drainvolume = 0;
prob.supplyvolume = 0;

space = [];
space.degree = [15 15];
space.coarsedegree = [8 8];
space.splitdim = [2 1 2 1];

time = [];
time.tol = 1e-4;
time.tspan = [0 5.258];
time.method = "NKS";

load initcond_pcl0.mat
time.initstate = finalstate;

model = blinkmulti(prob,space,time);

%%
model = solve(model);

%%
t = times(model,100);
v = volume(model,t);
plot(t,v/v(1)-1)

