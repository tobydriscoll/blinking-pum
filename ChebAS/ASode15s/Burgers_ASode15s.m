domain = [0 1;0 1];
split_flag = [false false];
degs = [10 10];
cdegs = [9 9];
tol = 1e-8;
R = 160;

f = @(x,y,t) 3/4 - 1./(4*(1+exp((R*(4*y-4*x-t)/32))));
g = @(x,y,t) 3/4 + 1./(4*(1+exp((R*(4*y-4*x-t)/32))));

cheb_struct.domain = domain;
cheb_struct.degs = degs;
cheb_struct.split_flag = split_flag;
cheb_struct.cdegs = cdegs;
cheb_struct.tol = tol;

Tree = ChebPatch(cheb_struct);
Tree = Tree.split(1);
Tree.split(2);

F = PUchebfun(Tree);
%set mass matricies
for i=1:length(F.leafArray)
   M{i} = BurgersMassMatrix(F.leafArray{i});
end

F.Setvalues(@(x,y)f(x,y,0));
u0 = F.Getvalues();
F.Setvalues(@(x,y)g(x,y,0));
v0 = F.Getvalues();

y0 = [u0;v0];


%set interpolation matrices
setInterpMatrices(F);

odetol = 1e-3;

tspan = [0 0.5];

opt = odeset('mass',M,'reltol',odetol,'abstol',odetol,'jacobian',@(t,y,approx)BurgersJacobian(t,y,approx,R),'BDF','on');
[t,U] = ASode15s(@(t,y,Approx) BurgersEvaluation(Approx,t,y,R),tspan,y0,F,opt);

%opt = odeset('mass',M{1},'reltol',odetol,'abstol',odetol,'jacobian',@(t,y)BurgersJacobian(t,y,Tree,R),'BDF','on');
%[t,U] = ode15s(@(t,y) BurgersEvaluation(Tree,t,y,R),tspan,y0,opt);
