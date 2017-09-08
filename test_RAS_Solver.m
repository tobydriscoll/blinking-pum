overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
split_flag = [1 1];
tol = 1e-12;
cheb_length  = 33;
dim = [33 33];
domain1 = [-1 overlap;-1 1];
domain2 = [-overlap 1;-1 1];

overlap_in = [-overlap overlap];

x = linspace(-1,1,100)';
[X,Y] = ndgrid(x);

children{1} = ChebPatch(domain1,deg_in,split_flag,tol);
children{2} = ChebPatch(domain2,deg_in,split_flag,tol);

force = @(x) ones(length(x),1);
border = @(x) zeros(length(x),1);

Tree = PUPatch(domain,overlap_in,2*cheb_length,children,1,[]);


sol1 = zeros(length(points1),1);
sol2 = zeros(length(points2),1);

sol = [sol1;sol2];

Tree.sample(sol);
%We will solve lap u = 1 with zero boundary conditions

newsol = RASSolve(Tree,domain,force,border,sol);

Tree.sample(newsol);

F = Tree.evalfGrid({x x},1,0);

surf(X,Y,F);
E = abs(F-TRUE);
E = max(E(:));
title(sprintf('The error is: %g',E));
