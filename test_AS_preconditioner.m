overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
cdeg_in = [3 3];
split_flag = [1 1];
tol = 1e-8;
gmres_tol = 5e-9;
maxit = 1200;

ep = 5e-3;
c = 1;
a=1;
b=0;

% Defaults for this blog post
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

%East West South North %BVP
%L = @(u,x,y,dx,dy,dxx,dyy) ep*dxx+ep*dyy+dx*(y-x^2)-dy*(2*x*y);
%B = {@(u,x,y,dx,dy,dxx,dyy) u, @(u,x,y,dx,dy,dxx,dyy) u, @(u,x,y,dx,dy,dxx,dyy) u,@(u,x,y,dx,dy,dxx,dyy) u};
%force = @(x) zeros(size(x,1),1);
%border = {@(x)ones(size(x,1),1),@(x)ones(size(x,1),1),@(x)ones(size(x,1),1),@(x)ones(size(x,1),1)};

% f1 = @(x,y,p) (-1).*16.^p.*p.*((-1).*((-1)+x).*x).^((-2)+p).*((-1).*((-1)+y).*y) ...
%   .^((-2)+p).*(((-1)+p).*((-1)+y).^2.*y.^2+(-2).*((-1)+2.*p).*x.*(( ...
%   -1)+y).^2.*y.^2+(-2).*x.^3.*((-1)+p.*(1+(-2).*y).^2+(-2).*((-1)+y) ...
%   .*y)+x.^4.*((-1)+p.*(1+(-2).*y).^2+(-2).*((-1)+y).*y)+x.^2.*((-1)+ ...
%   (-2).*((-1)+y).*y.*(1+((-1)+y).*y)+p.*(1+2.*((-1)+y).*y).^2));

%fb = @(x,y) (1-exp(-(1-x)/ep)).*(1-exp(-(1-y)/ep)).*cos(pi(x+y));

%f2 = @(x,y) exp(-2/ep)/ep*(-(exp(1/ep)-exp(y/ep)).*(-2*exp(1/ep)*ep^2*pi^2+exp(x/ep).*(1+2*ep^2*pi^2)).*cos(pi*(x+y))-(3*exp(2/ep)-exp((1+x)/ep)-exp((1+y)/ep)-exp((x+y)/ep)).*sin(pi*(x+y))*ep*pi);

% f2 = @(x,y) (-1).*exp(1).^((-1).*ep.^(-1).*(1+(-1).*y)).*(1+(-1).*exp(1).^(( ...
%  -1).*ep.^(-1).*(1+(-1).*x))).*ep.^(-1).*cos(pi.*(x+y))+(-1).*(1+( ...
%   -1).*exp(1).^((-1).*ep.^(-1).*(1+(-1).*x))).*(1+(-1).*exp(1).^(( ...
%   -1).*ep.^(-1).*(1+(-1).*y))).*pi.*sin(pi.*(x+y))+2.*((-1).*exp(1) ...
%   .^((-1).*ep.^(-1).*(1+(-1).*x)).*(1+(-1).*exp(1).^((-1).*ep.^(-1) ...
%   .*(1+(-1).*y))).*ep.^(-1).*cos(pi.*(x+y))+(-1).*(1+(-1).*exp(1).^( ...
%   (-1).*ep.^(-1).*(1+(-1).*x))).*(1+(-1).*exp(1).^((-1).*ep.^(-1).*( ...
%   1+(-1).*y))).*pi.*sin(pi.*(x+y)))+ep.*(exp(1).^((-1).*ep.^(-1).*( ...
%   1+(-1).*y)).*(1+(-1).*exp(1).^((-1).*ep.^(-1).*(1+(-1).*x))).* ...
%  ep.^(-2).*cos(pi.*(x+y))+exp(1).^((-1).*ep.^(-1).*(1+(-1).*x)).*( ...
%   1+(-1).*exp(1).^((-1).*ep.^(-1).*(1+(-1).*y))).*ep.^(-2).*cos(pi.* ...
%   (x+y))+2.*(1+(-1).*exp(1).^((-1).*ep.^(-1).*(1+(-1).*x))).*(1+(-1) ...
%   .*exp(1).^((-1).*ep.^(-1).*(1+(-1).*y))).*pi.^2.*cos(pi.*(x+y))+( ...
%   -2).*exp(1).^((-1).*ep.^(-1).*(1+(-1).*y)).*(1+(-1).*exp(1).^((-1) ...
%   .*ep.^(-1).*(1+(-1).*x))).*ep.^(-1).*pi.*sin(pi.*(x+y))+(-2).*exp( ...
%   1).^((-1).*ep.^(-1).*(1+(-1).*x)).*(1+(-1).*exp(1).^((-1).*ep.^( ...
%   -1).*(1+(-1).*y))).*ep.^(-1).*pi.*sin(pi.*(x+y)));
%   
% f3 = @(x,y) (1-exp((x-1)/ep)).*(1-exp((y-1)/ep)).*cos(pi*(x+y));

xc = -1;
yc = -1;
alpha = 10;
r0 = 1;
% 
f2 = @(x,y) ((x+(-1).*xc).^2+(y+(-1).*yc).^2).^(-1/2).*(alpha+alpha.^3.*( ...
r0.^2+(-1).*x.^2+2.*x.*xc+(-1).*xc.^2+(-1).*y.^2+2.*y.*yc+(-1).* ...
 yc.^2)).*((-1)+alpha.^2.*((-1).*r0.^2+(-1).*x.^2+2.*x.*xc+(-1).* ...
  xc.^2+(-1).*y.^2+2.*r0.*(x.^2+(-2).*x.*xc+xc.^2+(y+(-1).*yc).^2) ...
  .^(1/2)+2.*y.*yc+(-1).*yc.^2)).^(-2);

f3 = @(x,y) atan(alpha.*((-1).*r0+((x+(-1).*xc).^2+(y+(-1).*yc).^2).^(1/2)));

%c = sqrt(2);
%b = 20;
%d = 0.8;

% f2 = @(x,y) 2.*b.*c.*(c.^2.*(1+b.^2.*d.^2)+(-2).*b.^2.*c.*d.*(x+y)+b.^2.*(x+y) ...
%   .^2).^(-2).*((-1).*c.^2.*(1+b.^2.*d.^2).*((-1)+x+y).*((-1).*x+(-1) ...
%   .*y+2.*x.*y)+(-1).*b.^2.*(x+y).*((-1).*((-1)+x).*x.^2+x.^2.*((-3)+ ...
%   2.*x).*y+((-1)+x).*((-1)+2.*x).*y.^2+((-1)+2.*x).*y.^3)+2.*b.^2.* ...
%   c.*d.*((-1).*((-1)+x).*x.^2+x.*(1+2.*((-2)+x).*x).*y+((-1)+x).*(( ...
%   -1)+3.*x).*y.^2+((-1)+2.*x).*y.^3))+2.*(x+(-1).*x.^2+y+(-1).*y.^2) ...
%   .*atan(b.*c.^(-1).*((-1).*c.*d+x+y));
% 
% f3 = @(x,y) zeros(size(x,1),1);

%f2 = @(x,y) exp(x+y).*(y./(0.1+x)-2);

%f3 = @(x,y) exp(x+y);

%f2 = @(x,y) exp(1).^(0.109956E2.*y).*x.*((-6)+3.*exp(1).^(3.*x.*y).*x+(-0.119903E3).*x.^2+0.109956E2.*x.^2.*sin(pi.*(x+y)));

%f3 = @(x,y) x.^3.*exp(3.5*pi*y);

% f2 = @(x,y) exp(1).^((-0.1E3).*(0.263689E0+(-0.1E1).*x+x.^2+(-0.234E0).*y+ ...
% y.^2)).*(x.^3.*(0.8E5+(-0.8E5).*y).*y+(0.198E3+(-0.198E3).*y).*y+ ...
%   x.^4.*y.*((-0.4E5)+0.4E5.*y)+x.*(0.448E2+0.885396E4.*y+( ...
% -0.179076E5).*y.^2+0.4936E5.*y.^3+(-0.4E5).*y.^4)+x.^2.*(( ...
% -0.448E2)+(-0.48854E5).*y+0.579076E5.*y.^2+(-0.4936E5).*y.^3+ ...
% 0.4E5.*y.^4));
% 
% f3 = @(x,y) exp(1).^((-100).*(((-0.5E0)+x).^2+((-0.117E0)+y).^2)).*((-1)+x).* ...
%   x.*((-1)+y).*y;

%f2 = @(x,y) zeros(size(x,1),1);

%f3 = @(x,y) (cosh(10*x)+cosh(10*y))/(2*cosh(10));

%L = @(u,x,y,dx,dy,dxx,dyy) -ep*(dxx+dyy)+2*dx+dy;

L = @(u,x,y,dx,dy,dxx,dyy) (dxx+dyy);
B = {@(u,x,y,dx,dy,dxx,dyy) u, @(u,x,y,dx,dy,dxx,dyy) u, @(u,x,y,dx,dy,dxx,dyy) u,@(u,x,y,dx,dy,dxx,dyy) u};
force = @(x) f2(x(:,1),x(:,2));
border = {@(x)f3(x(:,1),x(:,2)),@(x)f3(x(:,1),x(:,2)),@(x)f3(x(:,1),x(:,2)),@(x)f3(x(:,1),x(:,2))};


cheb_struct.domain = domain;
cheb_struct.deg_in = deg_in;
cheb_struct.split_flag = split_flag;
cheb_struct.cdeg_in = cdeg_in;
cheb_struct.tol = tol;

Tree = ChebPatch(cheb_struct);

is_refined = false;

tic;

gmres_it = [];

ERRS = [];

NUMPTS = [];

while ~is_refined
    
    NUMPTS = [NUMPTS length(Tree)];
        
    if Tree.is_leaf
        
        [rhs] = setLinOps(Tree,L,B,force,border);
        
        sol = Tree.linOp\rhs;
        
        Max = Tree.sample(sol);
        
        
        Tree = Tree.splitleaf(Max,true);
    else
        
        [rhs] = setLinOps(Tree,L,B,force,border);
        Mat = CoarseASMat( Tree,L,B );
        
        A = @(sol) LaplacianForward(Tree,sol);
%        M = @(rhs) ASPreconditioner(Tree,rhs);
        M = @(rhs) CoarseCorrection(rhs,Tree,Mat);
        
        [sol,~,~,~,rvec] = gmres(A,rhs,[],gmres_tol,maxit,M,[],Tree.Getvalues);
        
        gmres_it = [gmres_it length(rvec)];
        
        length(rvec)
        
        Max = Tree.sample(sol);
        
        Tree.PUsplit(Max,true);
    end
    
    Tree.clean();
    
    x = linspace(domain(1,1),domain(1,2),100)';
    y = linspace(domain(2,1),domain(2,2),100)';
    [X,Y] = ndgrid(x,y);
    clf();
    
    G = Tree.evalfGrid({x y});
    subplot(1,2,1);
    surf(X,Y,G);
    xlabel('x');
    ylabel('y');
    subplot(1,2,2);
    Tree.plotdomain;
    pause(0.1);
    
    E = f3(X,Y) - G;
    ERRS = [ERRS norm(E(:),inf)];

    is_refined = Tree.is_refined;
end
toc



