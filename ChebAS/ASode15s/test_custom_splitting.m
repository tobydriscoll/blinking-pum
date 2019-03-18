domain = [-1 1;-1 1];
cheb_struct.domain = domain;
cheb_struct.degs = [20 20];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

odetol = 1e-3;
tspan = [0 0.3];

pctClosed = 0.75;

pA = 1e-4;
pS = 1e-3;
he = 2;

BoundaryH = 13;

%Test with 4 patches
Tree = ChebPatch(cheb_struct);

overlap = Tree.overlap;

Tree = Tree.split(1);
Tree.split(2);
 
Tree.split(1,false,0.08);
Tree.split(2,false,0.08);

Tree.reset();

H = PUchebfun(Tree);