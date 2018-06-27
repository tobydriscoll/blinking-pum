domain = [0 1;0 1];
split_flag = [false false];
degs = [32 32];
cdegs = [9 9];
tol = 1e-8;
R = 80;

f = @(x,y,t) 3/4 - 1./(4*(1+exp((R*(4*y-4*x-t)/32))));
g = @(x,y,t) 3/4 + 1./(4*(1+exp((R*(4*y-4*x-t)/32))));

cheb_struct.domain = domain;
cheb_struct.degs = degs;
cheb_struct.split_flag = split_flag;
cheb_struct.cdegs = cdegs;
cheb_struct.tol = tol;

Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);
%Tree.split(2);

F = PUchebfun(Tree);
%set mass matricies
for i=1:length(F.leafArray)
   M{i} = BurgersMassMatrix(F.leafArray{i});
end

%get intial condition
F.Setvalues(@(x,y)f(x,y,0));
u0 = Tree.Getvalues();
F.Setvalues(@(x,y)g(x,y,0));
v0 = Tree.Getvalues();

y0 = [u0;v0];

%set interpolation matrices
setInterpMatrices(F);

odetol = 1e-3;



tspan = [0 0.5];

opt = odeset('mass',M,'reltol',odetol,'abstol',odetol,'jacobian',@(t,y,approx)BurgersJacobian(t,y,approx,R));
%opt = odeset('mass',M{i},'reltol',odetol,'abstol',odetol);
%[t,U] = ode15s(@(t,y) BurgersEvaluation(Tree,t,y,R),tspan,y0,opt);

tspan = [0 0.5];

%[t,U] = ASode15s(@(Approx,t,y) BurgersEvaluation(Approx,t,y,R),tspan,y0,F,2,opt);

RHS = [];
for i=1:length(F.leafArray)
   ind = (i-1)*2*length(F.leafArray{i})+(1:(2*length(F.leafArray{i})));
   RHS = [RHS;M{i}*y0(ind)];
end

dh = 0.007;
z = ParPreconditionedNewtonForwardTime(0+dh,y0,RHS,F,@(Approx,t,y) BurgersEvaluation(Approx,t,y,R),2,dh,@(t,y,approx)BurgersJacobian(t,y,approx,R),M);
