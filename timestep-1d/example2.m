% Demo of SASPIN using Backward Euler for viscous Burgers

% includes coarse correction (2-level method)

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
clf, plot(x,u0,'ko')


%%
normres = []; normstep = [];  numgm = [];
un = u0;
u = un;
% time stepping loop
for ts = 1:ceil(T/dt)
    
    % solve for the new value using plain Newton
    for k = 1:5
        
        [y,yjac] = coarse_correction(u,un,dt,BC,dom,coarse);
        [z,zjac] = local_corrections(u+y,un,dt,BC,dom);
        
        normres(k,ts) = norm(y+z);
        if normres(k,ts) < 1e-12, break, end
        
        % find overall Newton step by GMRES
        yjac('clear'); zjac('clear');
        tol = min(0.1,1e-13*norm(u)/normres(k));
        jac = @(v) yjac(v) + zjac( v+yjac(v) );
        [s,~,~,~,gmhist] = gmres(jac,-y-z,[],tol);
        normstep(k,ts) = norm(s);  numgm(k,ts) = length(gmhist) - 1;

        % update
        u = u+s;
    end
    un = u;
    plot(x,un,'o'), ylim(sort(BC)),  drawnow
end

%%
% finite difference Jacobians, for sanity check
normres = []; normstep = [];  numgm = [];
un = u0;
u = un;
% time stepping loop
for ts = 1:1
    
    % solve for the new value using plain Newton
    for k = 1:5
        
        g = @(u) coarse_correction(u,un,dt,BC,dom,coarse);
        f = @(u) g(u) + local_corrections(u+g(u),un,dt,BC,dom);
        r = f(u);
        J = fdjac(f,u,r);
        
        normres(k,ts) = norm(r);
        if normres(k) < 1e-12, break, end
        
        % find overall Newton step by GMRES
        s = -(J\r);
        normstep(k,ts) = norm(s); 

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