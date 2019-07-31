prob = [];
prob.percentclosed = 0;
prob.A = 2e-6;
prob.S = 6.92e-5;
prob.h_slideover = 2;
prob.h_boundary = 13;
prob.drainvolume = 0;
prob.supplyvolume = 0;

space = [];
space.degree = [16 16];
space.coarsedegree = [8 8];
space.splitdim = [2 1 2 1 2];

time = [];
time.tol = 1e-4;
time.tspan = [0 5.258];
time.method = paropt.solver;
time.use_parallel = true;

%load initcond_pcl0.mat
time.initstate = []; %finalstate;

model = blinkmulti(prob,space,time);



