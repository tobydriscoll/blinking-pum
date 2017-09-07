overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
split_flag = [1 1];
tol = 1e-6;
maxit = 50;
cheb_length  = 33;
dim = [33 33];
domain1 = [-1 overlap;-1 1];
domain2 = [-overlap 1;-1 1];
border = @(x) exp(x(:,1)).*sin(x(:,2));

overlap_in = [-overlap overlap];

children{1} = ChebPatch(domain1,deg_in,split_flag,tol);
children{2} = ChebPatch(domain2,deg_in,split_flag,tol);

Tree = PUPatch(domain,overlap_in,2*cheb_length,children,1,[]);

points1 = Tree.children{1}.points;
points2 = Tree.children{2}.points;

sol1 = zeros(length(points1),1);
sol2 = zeros(length(points2),1);

out_border1 = FindBoundaryIndex2D(dim,Tree.children{1}.domain(),domain);
sol1(out_border1) = border(points1(out_border1,:));

out_border2 = FindBoundaryIndex2D(dim,Tree.children{2}.domain(),domain);
sol2(out_border2) =  border(points2(out_border2,:));

b = [sol1;sol2];

LaplacianForward(Tree,dim,domain,[sol1;sol2]);
RASPreconditioner(Tree,dim,domain,cheb_length,[sol1;sol2]);

A = @(sol) LaplacianForward(Tree,dim,domain,sol);
M = @(sol) RASPreconditioner(Tree,dim,domain,cheb_length,sol);

x1 = gmres(A,b,[],tol,maxit,M);

x = linspace(-1,1,100)';
[X,Y] = ndgrid(x);

Tree.sample(x1);

F = Tree.evalfGrid({x x},1,0);

surf(X,Y,F);