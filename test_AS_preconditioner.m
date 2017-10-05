overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
split_flag = [1 1];
tol = 1e-6;
maxit = 1200;
cheb_length  = 33;
dim = [33 33];

force = @(x) ones(length(x),1);

border = @(x) zeros(length(x),1);

Tree = ChebPatch(domain,deg_in,split_flag,tol);

%Tree = Tree.split(0.1/4,1);
%Tree.split(0.1/4);

%Tree.split(0.1/2);
%Tree.split(0.1/2);

%Tree.split(0.1);

Tree = Tree.split(0.1,1);
Tree.split(0.1);
% 
Tree.split(0.1);
Tree.split(0.1);
% 
Tree.split(0.1);

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
end


LaplacianForward(Tree,domain,zeros(Tree.length(),1));

A = @(sol) LaplacianForward(Tree,domain,sol);


M = @(rhs) ASPreconditioner(Tree,domain,rhs);
%M = @(rhs) CourseCorrection(rhs,Tree,domain);

G = @(sol) M(A(sol));
MM = @(rhs) M(A(M(rhs)));

tic
sol = gmres(A,rhs,[],tol,maxit,M);
toc

x = linspace(-1,1,100)';
[X,Y] = ndgrid(x);

Tree.sample(sol);

Fxx = Tree.evalfGrid({x x},1,2);
Fyy = Tree.evalfGrid({x x},2,2);

RESIDUAL = Fxx(:,:,3)+Fyy(:,:,3) - ones(100,100);
RESIDUAL([1 end],:)=Fxx([1 end],:,1);
RESIDUAL(:,[1 end])=Fxx(:,[1 end],1);
surf(X,Y,RESIDUAL);
figure();
surf(X,Y,Fxx(:,:,1));

