domain = [0 1;0 1];
cheb_struct.domain = domain;
cheb_struct.degs = [20 20];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

R = 10000;

odetol = 1e-4;
tspan = [0 1];

%Test with 2 patches
 Tree = ChebPatch(cheb_struct);
 Tree = Tree.split(1);
 Tree.split(2);

domain = [-1 1;-1 1];
cheb_struct.domain = domain;
cheb_struct.degs = [20 20];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

odetol = 1e-4;
tspan = [0 1];

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
%Tree.split(2);
Tree.split(1,false,0.1);
Tree.split(2,false,0.1); 
Tree.clean();

U_t = PUchebfun(Tree);
U_t.sample(@(x,y) zeros(size(x)));

V_t = U_t.copy();

setInterpMatrices(U_t,false);
setInterpMatrices(V_t,false);

[BurgerOp,M,y0] = SetBurgers(U_t,V_t,R,false);

opt = odeset('mass',M,'reltol',odetol,'abstol',odetol,'InitialStep',1);
[t,U] = ASode15s(false,BurgerOp,tspan,y0,{U_t,V_t},1,opt);