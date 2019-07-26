addpath ..
addpath ../PUChebfun

prob = [];
prob.percentclosed = 0;
prob.A = 2e-3;
prob.S = 1e-4;
prob.h_slideover = 2;
prob.h_boundary = 13;

space = [];
space.degree = [15 15];
space.coarsedegree = [8 8];
space.splitdim = [2 1 2];

time = [];
time.tol = 1e-4;
time.tspan = [0 5.258];
time.method = "NKS";

load initcond_pcl0.mat
time.initstate = finalstate;

model = blinkmulti(prob,space,time);

model = solve(model);

