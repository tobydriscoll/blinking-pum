overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
split_flag = [1 1];
tol = 1e-8;
maxit = 500;
cheb_length  = 33;
dim = [33 33];
domain1 = [-1 overlap;-1 1];
domain2 = [-overlap 1;-1 1];

force = @(x) ones(length(x),1);

border = @(x) zeros(length(x),1);

overlap_in = [-overlap overlap];

children{1} = ChebPatch(domain1,deg_in,split_flag,tol);
children{2} = ChebPatch(domain2,deg_in,split_flag,tol);

Tree = PUPatch(domain,overlap_in,2*cheb_length,children,1,[]);

points1 = Tree.children{1}.points;
points2 = Tree.children{2}.points;

sol1 = force(points1);
sol2 = force(points2);

out_border1 = FindBoundaryIndex2D(dim,Tree.children{1}.domain(),domain);
sol1(out_border1) = border(points1(out_border1,:));

out_border2 = FindBoundaryIndex2D(dim,Tree.children{2}.domain(),domain);
sol2(out_border2) =  border(points2(out_border2,:));

b = [sol1;sol2];

LaplacianForward(Tree,dim,domain,[sol1;sol2]);
F1 = border(points1);
F2 = border(points2);

A = @(sol) LaplacianForward(Tree,dim,domain,sol);


M = @(sol) RASStep(Tree,domain,force,border,sol);
G = @(sol) M(A(sol));
%sol = gmres(A,b,[],tol,maxit);
%sol = gmres(G,M(b),[],tol,maxit);
sol = gmres(A,b,[],tol,maxit,M);

x = linspace(-1,1,100)';
[X,Y] = ndgrid(x);

Tree.sample(sol);

F = Tree.evalfGrid({x x},1,0);

surf(X,Y,F);

R = A(sol)-b;

max(abs(R(:)))