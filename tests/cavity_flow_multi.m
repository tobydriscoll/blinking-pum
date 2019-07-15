addpath ../PUChebfun
addpath ../ChebAS
addpath ../ChebAS/examples

domain = [0 1;0 1];

F = PUchebfun(@(x,y)exp(-y.^20./(1-y.^20)),[0 1;0 1],'Degree',[33 33],'CoarseDegree',[9 9],'tol',1e-5,'SplitAll',true); F.reset()
F.sample(@(x,y) zeros(size(x)));

setInterpMatrices(F,true);

steep = 0.1;

%F.Setvalues(@(x,y)SideBumpFunc(y,[0 1],steep)); %set values
%F.sample(@(x,y)SideBumpFunc(y,[0 1],steep)); %set coeffs
bump = @(x,y) y.*(1 - SideBumpFunc(x,[0 1],steep) - SideBumpFunc(1-x,[0 1],steep));
F.Setvalues(bump); %set values
F.sample(bump); %set coeffs

u = F.Getvalues();
Fy = diff(F,2,1);
w = -(Fy.Getvalues());
v = zeros(length(F),1);
init = [u;v;w];

opt = [];
opt.numcomp = 3;  % # of solution components
opt.reltol = 1e-10;
opt.abstol = 1e-11;
opt.inparallel = true;
opt.numwork = 2;   % number of workers
opt.coarsetol = 1e-4;
opt.griddepth = 2;
opt.twolevel = true;
%opt.solver = @NKSsolver;
opt.solver = @SNKsolver;

Re = [100 500 750 1000];

for kk = 1:length(Re) 
    f = @(u,leaf) CavityFlow(Re(kk),u,leaf,steep);
    Jac = @(u,leaf) CavityFlowJacobian(Re(kk),u,leaf);
    
    
    if opt.inparallel
        pp = gcp('nocreate');
        if isempty(pp) || pp.NumWorkers ~=opt.numwork
            delete(pp);
            pp = parpool(opt.numwork);
        end
        
        bstart = pp.ticBytes;
        tstart = tic;
        [ sol,stats ] = opt.solver(f,Jac,init,F,opt);
        stats.elapsedtime = toc(tstart);
        stats.tocbytes = pp.tocBytes(bstart);
    else
        tstart = tic;
        [ sol,stats ] = opt.solver(f,Jac,init,F,opt);
        stats.elapsedtime = toc(tstart);
    end
    
    results{kk} = {opt Re(kk) stats F sol};
    save cavity_flow_1 results
    
    init = sol;
    
end