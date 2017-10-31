overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
cdeg_in = [3 3];
split_flag = [1 1];
tol = 1e-7;
gmres_tol = 1e-7;
maxit = 1200;
cheb_length  = 65;
dim = [65 65];

c = 1e-2;
a=1;
b=0;

L = @(u,x,y,dx,dy,dxx,dyy) (dxx+dyy);

%East West South North
B = {@(u,x,y,dx,dy,dxx,dyy) u, @(u,x,y,dx,dy,dxx,dyy) u, @(u,x,y,dx,dy,dxx,dyy) u, @(u,x,y,dx,dy,dxx,dyy) u};

%force = @(x) 4*c./(c+x(:,1).^2+x(:,2).^2).^2;
%border = @(x) log((x(:,1).^2+x(:,2).^2)/c+1);

force = @(x) ones(length(x),1);
border = @(x) zeros(length(x),1);

Tree = ChebPatch(domain,domain,domain,deg_in,split_flag,tol,cdeg_in);

is_refined = false;

tic;
while ~is_refined
    
    if Tree.is_leaf
        
        [rhs] = setLinOps(Tree,L,B,force,border);
        
        sol = Tree.linOp\rhs;
        
        Tree.sample(sol);
        
        
        Tree = Tree.splitleaf(true);
    else
        
        [rhs] = setLinOps(Tree,L,B,force,border);
        Mat = CoarseASMat( Tree,L,B );
        
        A = @(sol) LaplacianForward(Tree,domain,sol);
        %M = @(rhs) ASPreconditioner(Tree,domain,rhs);
        M = @(rhs) CoarseCorrection(rhs,Tree,domain,Mat);
        
        [sol,~,~,~,rvec] = gmres(A,rhs,[],gmres_tol,maxit,M,[],Tree.Getvalues);
        
        Tree.sample(sol);
        
        Tree.PUsplit(true);
    end
    
    x = linspace(-1,1,100)';
    [X,Y] = ndgrid(x);
    clf();
    
    G = Tree.evalfGrid({x x},1,0);
    subplot(1,2,1);
    surf(X,Y,G);
    subplot(1,2,2);
    Tree.plotdomain;
    pause(0.1);
    
    is_refined = Tree.is_refined;
end
toc



