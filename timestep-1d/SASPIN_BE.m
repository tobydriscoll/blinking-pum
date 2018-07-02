% Demo of SASPIN using Backward Euler for viscous Burgers

% define the problem
BC = [1 -1.4];
ODE = struct;
ODE.fun = @(u,D1,D2) 0.02*(D2*u) - u.*(D1*u);
ODE.jac = @(u,D1,D2) 0.02*D2 - u.*D1 - diag(D1*u);

% define the subdomains
dom = struct;
%dom(1).xlim = [-1 0.06];  dom(1).n = 50;
%dom(2).xlim = [-0.05 1];  dom(2).n = 47;
dom(1).xlim = [-1 -0.48];  dom(1).n = 60;
dom(2).xlim = [-0.52 0.6];  dom(2).n = 100;
dom(3).xlim = [0.54 1];  dom(3).n = 50;
n = cat(1,dom.n);

% time stepping parameters
T = 1;  dt = .02;

[data,DAEfun,jacfun] = setup_BE(ODE,BC,dt,dom);

% initial condition
x = cat(1,data.x);
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
        [z,jacfun] = corrections(u,un,dt,BC,dom,data);
        normres(k,ts) = norm(z);
        
        if normres(k) < 1e-15, break, end
        
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