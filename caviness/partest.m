addpath ../../chebfun
addpath ../PUChebfun
addpath ..

nw_ = [1 2 4 6 8 12 16 20];

p = []; 
solver = string([]);
elapsed = [];
numouter = [];
numinner = [];
bytesout = [];
bytesin = [];
restime = [];
jactime = [];

for nw = nw_
    paropt.numwork = nw;
    pp = gcp('nocreate');
    if isempty(pp) || pp.NumWorkers ~=paropt.numwork
        delete(pp);
        pp = parpool(paropt.numwork);
    end
    
    for method = ["SNK" "NKS"]
      paropt.solver = method;
	  setup_model
	  
      bstart = pp.ticBytes;
      tstart = tic;	  
	  model = solve(model);
      stats.elapsedtime = toc(tstart);
      stats.tocbytes = pp.tocBytes(bstart);

	  p(end+1,1) = nw;
	  solver(end+1,1) = method;
	  elapsed(end+1,1) = stats.elapsedtime;
	  %restime(end+1,1) = stats.restime;
	  %jactime(end+1,1) = stats.jactime;
	  %numouter(end+1,1) = length(stats.normres)-1;
	  %numinner(end+1,1) = sum(stats.numgmres);
	  bytesin(end+1,1) = mean(stats.tocbytes(:,1),1);
	  bytesout(end+1,1) = mean(stats.tocbytes(:,2),1);
    end

end

%results = table(p,solver,time,numouter,numinner,bytesin,bytesout)
results = table(p,solver,elapsed,bytesin,bytesout)
save partest_1_parfor results prob space time paropt

