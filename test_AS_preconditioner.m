overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
split_flag = [1 1];
tol = 1e-6;
maxit = 1200;
cheb_length  = 33;
dim = [33 33];

domain1 = [-1 overlap;-1 1];
domain2 = [-overlap 1;-1 1];

domain11 = [-1 overlap;-1 overlap];
domain12 = [-1 overlap;-overlap 1];

domain21 = [-overlap 1;-1 overlap];
domain22 = [-overlap 1;-overlap 1];

force = @(x) ones(length(x),1);

border = @(x) zeros(length(x),1);

overlap_in = [-overlap overlap];

childrenl{1} = ChebPatch(domain11,deg_in,split_flag,tol);
childrenl{2} = ChebPatch(domain12,deg_in,split_flag,tol);

childrenr{1} = ChebPatch(domain21,deg_in,split_flag,tol);
childrenr{2} = ChebPatch(domain22,deg_in,split_flag,tol);

children{1} = PUPatch(domain1,overlap_in,2*cheb_length,childrenl,2,[]);
children{2} = PUPatch(domain2,overlap_in,2*cheb_length,childrenr,2,[]);

Tree = PUPatch(domain,overlap_in,4*cheb_length,children,1,[]);

%children{1} = ChebPatch(domain1,deg_in,split_flag,tol);
%children{2} = ChebPatch(domain2,deg_in,split_flag,tol);

%Tree = PUPatch(domain,overlap_in,2*cheb_length,children,1,[]);

LEAVES = Tree.collectLeaves({});

rhs = [];

for k=1:length(LEAVES)
    points = LEAVES{k}.points;
    sol = force(points);
    dim = LEAVES{k}.degs;
        
    Dxx = kron(eye(dim(2)),diffmat(dim(1),2));
    Dyy = kron(diffmat(dim(2),2),eye(dim(1)));
    
    scalex = 2/diff(LEAVES{k}.domain(1,:));
    scaley = 2/diff(LEAVES{k}.domain(2,:));
    
    lap = scalex^2*Dxx+scaley^2*Dyy;
    
    [out_border,in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
    interior = ~out_border & ~ in_border;
    sol(out_border) = border(points(out_border,:));
    sol(in_border) = 0;
    rhs = [rhs;sol];
end


LaplacianForward(Tree,domain,zeros(Tree.length(),1));

A = @(sol) LaplacianForward(Tree,domain,sol);


M = @(rhs) ASPreconditioner(Tree,domain,rhs);
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
RESIDUAL([1 end],:)=0;
RESIDUAL(:,[1 end])=0;
surf(X,Y,RESIDUAL);
figure();
surf(X,Y,Fxx(:,:,1));