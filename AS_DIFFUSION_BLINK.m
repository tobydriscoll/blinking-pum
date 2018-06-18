overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
cdeg_in = [3 3];
split_flag = [1 1];
tol = 1e-4;
gmres_tol = 5e-5;
maxit = 1200;
T = 5;
dt = 0.001;
theta = 0.5;
t0 = 0.01;
x0 = 0.1;
y0 = 0.2;
pc = 0.8;


alpha = 1.6;  gamma = 7;

freq = 2*16*pi;

lambda = @(t) 1 + pc*(-1+tanh(4*cos(freq*t)));  % upper eyelid position, in (-1,1]
dlambda_dt = @(t) freq*4*pc*sin(freq*t).*sech(4*cos(freq*t)).^2;


dhx = @(dx,x) ((x.^2-alpha^2).^2)./((x.^2+alpha^2)*gamma)*dx;
dhy = @(dy,t) 2./(lambda(t)+1)*dy;

y_s = @(y,t) 1/2*(y+1).*(lambda(t)+1)-1;
x_s = @(x) gamma*x./(alpha^2-x.^2);

y_e = @(x,y) sin(y)./(cos(y)+cosh(x));
x_e = @(x,y) sinh(x)./(cos(y)+cosh(x));


%EAST WEST SOUTH NORTH
%B = {@(h,x,y,dx,dy,dxx,dyy,t) dx, @(h,x,y,dx,dy,dxx,dyy,t) dx, @(h,x,y,dx,dy,dxx,dyy,t) dy,@(h,x,y,dx,dy,dxx,dyy,t) dy};
B = {@(u,x,y,dx,dy,dxx,dyy,t) u, @(u,x,y,dx,dy,dxx,dyy,t) u, @(u,x,y,dx,dy,dxx,dyy,t) u,@(u,x,y,dx,dy,dxx,dyy,t) u};

%bf = @(x,y,t) zeros(size(x,1),1);


kernal = @(x,y,t) exp( -(((x-x0).^2+(y-y0).^2)/(4*(t+t0)))) / (4*pi*(t+t0));
bf = @(x,y,t) kernal(x_e(x_s(x),y_s(y,t)),y_e(x_s(x),y_s(y,t)),t);
border = {bf,bf,bf,bf};
init = @(x,y) bf(x,y,0);
exact = @(x,y,t) bf(x,y,t);

cheb_struct.domain = domain;
cheb_struct.deg_in = deg_in;
cheb_struct.split_flag = split_flag;
cheb_struct.cdeg_in = cdeg_in;
cheb_struct.tol = tol;

Tree = ChebPatch(cheb_struct);
% Tree = Tree.split(1);
% Tree.split(2);

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

%plot(F);

x = linspace(-1,1,100)';
[X,Y] = ndgrid(linspace(-1,1,100));

X_s = x_s(X); Y_s = y_s(Y,CT);

Z = tanh((X_s+1i*Y_s)/2);

X_e = sinh(X_s)./(cos(Y_s)+cosh(X_s));  Y_e = sin(Y_s)./(cos(Y_s)+cosh(X_s));




while CT<T
    if length(F.leafArray)>1
        rhs = setLinOpsThetaBlink(F,B,border,theta,dt,CT,lambda,dlambda_dt,alpha,gamma);
        Mat = CoarseASMatThetaBlink(F,B,theta,dt,CT,lambda,dlambda_dt,alpha,gamma);
        
        A = @(sol) ParSchwarzForward(F,sol);
        M = @(rhs) CoarseCorrection(F,rhs,Mat);
        
        tic, [sol,~,~,~,rvec] = gmres(A,rhs,[],gmres_tol,maxit,M,[],F.ChebRoot.Getvalues); toc;
        num_it = length(rvec);
        
    else
        rhs = setLinOpsThetaBlink(F,B,border,theta,dt,CT,lambda,dlambda_dt,alpha,gamma);
        
        sol = F.ChebRoot.linOp\rhs;
        
        num_it = 0;
    end
    F.ChebRoot.sample(sol);
    
    
    CT = CT+dt;
    
    
    
    
    X_s = x_s(X); Y_s = y_s(Y,CT);
    
    Z = tanh((X_s+1i*Y_s)/2);
    
    X_e = sinh(X_s)./(cos(Y_s)+cosh(X_s));  Y_e = sin(Y_s)./(cos(Y_s)+cosh(X_s));
    
    F_exact = exact(X,Y,CT);
    
    F_eval = F.evalfGrid({x x});
    Error = F_exact-F_eval;
    
    Error = norm(Error(:),inf);
    DT = delaunay(X_e(:),Y_e(:));
    
    
    
    subplot(1,2,1);
pcolor(X_e,Y_e,F_eval);
    title(sprintf('time %g iterations %d relative error %2.2e',CT,num_it,Error));
xlim([-1 1]);
ylim([-1 1]);
view([360 90]);
    caxis([0 8])
shading flat;
subplot(1,2,2);
pcolor(X_e,Y_e,F_exact);
title('exact');
xlim([-1 1]);
ylim([-1 1]);
view([360 90]);
    caxis([0 8])
shading flat;

pause(0.001);
end




