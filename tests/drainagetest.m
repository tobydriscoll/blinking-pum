addpath ..
addpath ../PUChebfun

prob = [];
prob.percentclosed = 0;
prob.A = 2e-6;
prob.S = 6.92e-5;
prob.h_slideover = 2;
prob.h_boundary = 13;
prob.drainvolume = 8;
prob.supplyvolume = 8;

space = [];
space.degree = [16 16];
space.coarsedegree = [8 8];
space.splitdim = [2 1 2 1 2];

time = [];
time.tol = 1e-4;
time.tspan = [0 5.258];
time.method = "SNK";
time.use_parallel = false;

load initcond_pcl0.mat
time.initstate = finalstate;

model = blinkmultilog(prob,space,time);

%%
model = solve(model);

% %%
% t = times(model,100);
% v = volume(model,t);
% plot(t,v/v(1)-1)
% xlabel('time')
% ylabel('relative volume change')

%%
animate(model,'drainage_8_pcl0.mp4',linspace(0,5.258,601))

%%
save drainage_8_pcl0 model