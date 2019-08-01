addpath ../../chebfun
addpath ../PUChebfun
addpath ..

nw = 16
    paropt.numwork = nw;
    pp = gcp('nocreate');
    if isempty(pp) || pp.NumWorkers ~=paropt.numwork
        delete(pp);
        pp = parpool(paropt.numwork)
    end
    
prob = [];
prob.percentclosed = 0;
prob.A = 2e-6;
prob.S = 6.92e-5;
prob.h_slideover = 2;
prob.h_boundary = 13;
prob.drainvolume = 8;
prob.supplyvolume = 8;

space = [];
space.degree = [20 20];
space.coarsedegree = [8 8];
space.splitdim = [2 1 2 1 2 1];

time = [];
time.tol = 1e-4;
time.tspan = [0 5.258];
time.method = "SNK";
time.use_parallel = true;

load initcond_pcl0.mat
time.initstate = finalstate;

model = blinkmultilog(prob,space,time);

	  
      bstart = pp.ticBytes
      tstart = tic	  
	  model = solve(model);
      elapsedtime = toc(tstart)
      tocbytes = pp.tocBytes(bstart)

save drainage_64_pcl0 model
