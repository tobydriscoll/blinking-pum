domain = [0 1;0 1];
cheb_struct.domain = domain;
cheb_struct.degs = [30 3];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

R = 800;

odetol = 1e-4;
tspan = [0 1];

%Test with 2 patches
%  Tree = ChebPatch(cheb_struct);
%  Tree = Tree.split(1);
%  Tree.split(2);

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
%Tree.split(2);
Tree.split(1,false,0.1);
Tree.split(2,false,0.1); 

U_t = PUchebfun(Tree);
U_t.sample(@(x,y) zeros(size(x)));

V_t = U_t.copy();

setInterpMatrices(U_t,false);
setInterpMatrices(V_t,false);

[BurgerOp,M,y0] = SetBurgers(U_t,V_t,R,false);


yp0 = Masstimes({U_t,V_t},M,ParLocalResidual(0,y0,1,{U_t,V_t},BurgerOp));


%[y,yp] = GetInitialSlope(M,y0,zeros(length(y0),1),0,{U_t,V_t},BurgerOp,1e-3);
opt = odeset('mass',M,'reltol',odetol,'abstol',odetol);
[t,U] = ASode15s(true,BurgerOp,tspan,y0,{U_t,V_t},opt);