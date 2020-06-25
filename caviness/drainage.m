addpath ../PUChebfun
addpath ..

load initcond_pcl0.mat
initstate = finalstate
%disp(min2(finalstate.H))
%ratio = initstate.dH ./ initstate.H
%H = initstate.H;
%dH = initstate.dH;
%f = @(x,y) dH(x,y)./H(x,y)
%ratio = chebfun2(f)

% nw = 16
% myCluster=parallel.cluster.Local;
% myCluster.NumWorkers=nw;
% %myCluster.NumThreads=2;
% myCluster.JobStorageLocation=getenv('TMPDIR');
% %setenv('OMP_THREAD_LIMIT', int2str(myCluster.NumThreads));
% pp = parpool(myCluster);

%    paropt.numwork = nw;
%    pp = gcp('nocreate');
%    if isempty(pp) || pp.NumWorkers ~=paropt.numwork
%        delete(pp);
%        pp = parpool(paropt.numwork)
%    end
    
prob = [];
prob.percentclosed = 0;
prob.A = 2e-6;
prob.S = 6.92e-5;
prob.h_slideover = 2;
prob.h_boundary = 13;
prob.drainvolume = 8;
prob.supplyvolume = 8;

space = [];
space.degree = [24 32];
space.coarsedegree = [8 8];
space.splitdim = [2 1 2 1];

time = [];
time.tol = 1e-4;
time.tspan = [0 0.005];%[0 5.258];
time.method = "SNK";
time.use_parallel = false;
time.initstate = initstate;

model = blinkmulti(prob,space,time)
disp('model created')
	  
%bstart = pp.ticBytes;
tstart = tic
model = solve(model);
elapsedtime = toc(tstart)
%tocbytes = pp.tocBytes(bstart)

%save drainage_8_nw16_pcl0 model
