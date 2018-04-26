overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
cdeg_in = [3 3];
split_flag = [1 1];
tol = 1e-8;
gmres_tol = 5e-9;
maxit = 1200;
T = 5;
dt = 0.1;
theta = 0.5+dt;

alpha = 1;
L = @(u,x,y,dx,dy,dxx,dyy,t) alpha*(dxx+dyy);

%EAST WEST SOUTH NORTH
%B = {@(u,x,y,dx,dy,dxx,dyy,t) -dx, @(u,x,y,dx,dy,dxx,dyy,t) dx, @(u,x,y,dx,dy,dxx,dyy,t) -dy,@(u,x,y,dx,dy,dxx,dyy,t) dy};
B = {@(u,x,y,dx,dy,dxx,dyy,t) u, @(u,x,y,dx,dy,dxx,dyy,t) u, @(u,x,y,dx,dy,dxx,dyy,t) u,@(u,x,y,dx,dy,dxx,dyy,t) u};
force = @(x,y) zeros(size(x,1),1);
time_dependent = false;

bf = @(x,y) zeros(size(x,1),1);
border = {bf,bf,bf,bf};
init = @(x,y) exp(-1./(1-x.^2)).*exp(-1./(1-y.^2));
%init = @(x,y) x.*y;

cheb_struct.domain = domain;
cheb_struct.deg_in = deg_in;
cheb_struct.split_flag = split_flag;
cheb_struct.cdeg_in = cdeg_in;
cheb_struct.tol = tol;

Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);
Tree.split(2);

for k=1:1
    Tree.split(1);
    Tree.split(2);
end

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

%rhs = setLinOpsTheta(F,L,B,border,theta,dt,CT);

%Going to assume we have more than one patch
while CT<T
    rhs = setLinOpsTheta(F,L,B,border,theta,dt,CT);
    Mat = CoarseASMatTheta(F,L,B,theta,dt,CT);
    
    A = @(sol) ParSchwarzForward(F,sol);
    M = @(rhs) CoarseCorrection(F,rhs,Mat);
    
    tic, [sol,~,~,~,rvec] = gmres(A,rhs,[],gmres_tol,maxit,M,[],F.ChebRoot.Getvalues); toc;
    length(rvec)
    F.ChebRoot.sample(sol);
    sum(F^2)
    
    plot(F);
    
    CT = CT+dt;
    
    title(sprintf('time %g',CT));
    
    pause(0.001);
end


% while CT<T
%     
%     setLinOps(F.ChebRoot,L,B,force,border);
%     
%     
%     tic,AS_P_STEP(F.ChebRoot,dt,theta);toc
%     
%     plot(F);
%     
%     CT = CT+dt;
%     
%     title(CT);
%     
%     pause(0.001);
% end



