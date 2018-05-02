overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
cdeg_in = [3 3];
split_flag = [1 1];
tol = 1e-8;
gmres_tol = 5e-9;
maxit = 1200;
T = 1;
dt = 0.0001;
theta = 0.5+dt;

alpha = 0.1;
L = @(u,x,y,dx,dy,dxx,dyy,t) alpha*(dxx+dyy);

%EAST WEST SOUTH NORTH
B = {@(u,x,y,dx,dy,dxx,dyy,t) dx, @(u,x,y,dx,dy,dxx,dyy,t) dx, @(u,x,y,dx,dy,dxx,dyy,t) dy,@(u,x,y,dx,dy,dxx,dyy,t) dy};
%B = {@(u,x,y,dx,dy,dxx,dyy,t) u, @(u,x,y,dx,dy,dxx,dyy,t) u, @(u,x,y,dx,dy,dxx,dyy,t) u,@(u,x,y,dx,dy,dxx,dyy,t) u};
force = @(x,y) zeros(size(x,1),1);
time_dependent = false;

bf = @(x,y) zeros(size(x,1),1);
border = {bf,bf,bf,bf};
init = @(x,y) cos(pi*x).*cos(pi*y)+cos(4*pi*x).*cos(4*pi*y);
exact = @(x,y,t) exp(-alpha*2*pi^2*t).*cos(pi*x).*cos(pi*y)+exp(-alpha*32*pi^2*t).*cos(4*pi*x).*cos(4*pi*y);


cheb_struct.domain = domain;
cheb_struct.deg_in = deg_in;
cheb_struct.split_flag = split_flag;
cheb_struct.cdeg_in = cdeg_in;
cheb_struct.tol = tol;

Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);
Tree.split(2);

% for k=1:1
%     Tree.split(1);
%     Tree.split(2);
% end

is_refined = false;

tic;

gmres_it = [];

ERRS = [];

NUMPTS = [];

CT = 0;

Tree.sample(init);

F = PUchebfun(Tree);

sum(F^2)

plot(F);
pause(0.001);

x = linspace(-1,1,100)';
[X,Y] = ndgrid(linspace(-1,1,100));

while CT<T
    rhs = setLinOpsTheta(F,L,B,border,theta,dt,CT);
    Mat = CoarseASMatTheta(F,L,B,theta,dt,CT);
    
    A = @(sol) ParSchwarzForward(F,sol);
    M = @(rhs) CoarseCorrection(F,rhs,Mat);
    
    tic, [sol,~,~,~,rvec] = gmres(A,rhs,[],gmres_tol,maxit,M,[],F.ChebRoot.Getvalues); toc;
    num_it = length(rvec);
    F.ChebRoot.sample(sol);
    
    plot(F);
    
    CT = CT+dt;
    
    F_exact = exact(X,Y,CT);
    Error = F_exact - F.evalfGrid({x x});
    
    Error = norm(Error(:),inf);
    
    title(sprintf('time %g iterations %d relative error %2.2e',CT,num_it,Error));
    
    pause(0.001);
end




