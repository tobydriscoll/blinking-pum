overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
cdeg_in = [3 3];
split_flag = [1 1];
tol = 1e-6;
maxit = 1200;
cheb_length  = 33;
dim = [33 33];

force = @(x) ones(length(x),1);

border = @(x) zeros(length(x),1);

Tree = ChebPatch(domain,domain,domain,deg_in,split_flag,tol,cdeg_in);

Tree = Tree.split(1);
Tree.split();

Tree.split();
Tree.split();
 
Tree.split();

LEAVES = Tree.collectLeaves({});

rhs = [];

for k=1:length(LEAVES)
    points = LEAVES{k}.points;
    sol = force(points);
    dim = LEAVES{k}.degs;
    
    [out_border,in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
    
    interior = ~out_border & ~ in_border;
    
    sol(out_border) = border(points(out_border,:));
    sol(in_border) = 0;
    rhs = [rhs;sol];
    
    Dxx = kron(eye(dim(2)),diffmat(dim(1),2,LEAVES{k}.domain(1,:)));
    Dyy = kron(diffmat(dim(2),2,LEAVES{k}.domain(2,:)),eye(dim(1)));
    
    E = eye(prod(dim));
    
    Lap = Dxx+Dyy;
    
    Lap(in_border,:) = E(in_border,:);
    Lap(out_border,:) = E(out_border,:);
    
    LEAVES{k}.linOp = Lap;
    
    cdim = LEAVES{k}.cdegs;
    [out_border,in_border] = FindBoundaryIndex2D(cdim,LEAVES{k}.domain(),domain);
    
    Dxx = kron(eye(cdim(2)),diffmat(cdim(1),2,LEAVES{k}.domain(1,:)));
    Dyy = kron(diffmat(cdim(2),2,LEAVES{k}.domain(2,:)),eye(cdim(1)));
    
    E = eye(prod(cdim));
    
    CLap = Dxx+Dyy;
    
    CLap(in_border,:) = E(in_border,:);
    CLap(out_border,:) = E(out_border,:);
    
    LEAVES{k}.ClinOp = CLap;
end


LaplacianForward(Tree,domain,zeros(Tree.length(),1));

A = @(sol) LaplacianForward(Tree,domain,sol);


%M = @(rhs) ASPreconditioner(Tree,domain,rhs);
M = @(rhs) CoarseCorrection(rhs,Tree,domain);

tic
sol = gmres(A,rhs,[],tol,maxit,M);
toc

x = linspace(-1,1,100)';
[X,Y] = ndgrid(x);

Tree.sample(sol);

G = Tree.evalfGrid({x x});

surf(X,Y,G);


