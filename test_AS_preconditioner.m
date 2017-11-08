overlap = 0.1;

domain = [-2 2;-2 2];

deg_in = [5 5];
cdeg_in = [3 3];
split_flag = [1 1];
tol = 1e-8;
gmres_tol = tol;
maxit = 1200;
cheb_length  = 65;
dim = [65 65];

ep = 5e-2;
c = 1;
a=1;
b=0;

%East West South North %BVP
L = @(u,x,y,dx,dy,dxx,dyy) ep*dxx+ep*dyy+dx*(y-x^2)-dy*(2*x);
B = {@(u,x,y,dx,dy,dxx,dyy) u, @(u,x,y,dx,dy,dxx,dyy) u, @(u,x,y,dx,dy,dxx,dyy) u,@(u,x,y,dx,dy,dxx,dyy) u};
force = @(x) zeros(size(x,1),1);
border = {@(x)ones(size(x,1),1),@(x)ones(size(x,1),1),@(x)ones(size(x,1),1),@(x)ones(size(x,1),1)};

% L = @(u,x,y,dx,dy,dxx,dyy) dxx+dyy;
% B = {@(u,x,y,dx,dy,dxx,dyy) u, @(u,x,y,dx,dy,dxx,dyy) u, @(u,x,y,dx,dy,dxx,dyy) u,@(u,x,y,dx,dy,dxx,dyy) u};
% force = @(x) ones(size(x,1),1);
% border = {@(x)zeros(size(x,1),1),@(x)zeros(size(x,1),1),@(x)zeros(size(x,1),1),@(x)zeros(size(x,1),1)};

Tree = ChebPatch(domain,domain,domain,deg_in,split_flag,tol,cdeg_in);

is_refined = false;

tic;
while ~is_refined
    
    if Tree.is_leaf
        
        [rhs] = setLinOps(Tree,L,B,force,border);
        
        sol = Tree.linOp\rhs;
        
        Max = Tree.sample(sol);
        
        
        Tree = Tree.splitleaf(Max,true);
    else
        
        [rhs] = setLinOps(Tree,L,B,force,border);
        Mat = CoarseASMat( Tree,L,B );
        
        A = @(sol) LaplacianForward(Tree,domain,sol);
        %M = @(rhs) ASPreconditioner(Tree,domain,rhs);
        M = @(rhs) CoarseCorrection(rhs,Tree,domain,Mat);
        
        [sol,~,~,~,rvec] = gmres(A,rhs,[],gmres_tol,maxit,M,[],Tree.Getvalues);
        
        Max = Tree.sample(sol);
        
        Tree.PUsplit(Max,true);
    end
    
    x = linspace(domain(1,1),domain(1,2),100)';
    y = linspace(domain(2,1),domain(2,2),100)';
    [X,Y] = ndgrid(x,y);
    clf();
    
    G = Tree.evalfGrid({x y});
    subplot(1,2,1);
    surf(X,Y,G);
    xlabel('x');
    ylabel('t');
    subplot(1,2,2);
    Tree.plotdomain;
    pause(0.1);
    
    is_refined = Tree.is_refined;
end
toc



