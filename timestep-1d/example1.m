% Demo of SASPIN using Backward Euler for viscous Burgers

% define the problem
BC = [1 -1.4];
ODE = struct;
ODE.fun = @(u,D1,D2) 0.02*(D2*u) - u.*(D1*u);
ODE.jac = @(u,D1,D2) 0.02*D2 - u.*D1 - diag(D1*u);

% define the subdomains
Nd = 20;
dom = struct;
dx = 2/Nd;
zone = -1 + dx*(0:Nd);
Npt = 1000;  dn = round(Npt./[1.1*Nd 0.9*Nd]);
for i = 1:Nd
    r = dx*(.01 + 0.05*rand(1));
    xlim = zone(i:i+1);
    if i > 1, xlim(1) = xlim(1)-r;  end
    if i < Nd, xlim(2) = xlim(2)+r; end
    dom(i).xlim = xlim;
    dom(i).n = randi(dn);
end

n = cat(1,dom.n);

% time stepping parameters
T = 2;  dt = .02;

[dom,coarse] = setup(ODE,dom);

% initial condition
x = cat(1,dom.x);
u0 = BC(1) + (x+1)/2*diff(BC);
un = u0;
clf, plot(x,u0,'ko')


%%
normres = []; normstep = [];  numgm = [];
u = un;
% time stepping loop
for ts = 1:ceil(T/dt)
    
    % solve for the new value using plain Newton
    for k = 1:5
        % evaluate the local corrections/solve local nonlinear problems
        [z,jacfun] = local_corrections(u,un,dt,BC,dom);
        normres(k,ts) = norm(z);
        
        if normres(k) < 1e-12, break, end
        
        % find overall Newton step by GMRES
        jacfun('clear')
        tol = min(0.1,1e-10*norm(u)/normres(k));
        [s,~,~,~,gmhist] = gmres(jacfun,-z,[],tol);
        normstep(k,ts) = norm(s);  numgm(k,ts) = length(gmhist) - 1;

        % update
        u = u+s;
    end
    un = u;
    plot(x,un,'o'), ylim(sort(BC)),  drawnow
end

%%
% check the Newton residuals, GMRES steps
semilogy(normres')
maxGMiter = max(numgm(:))