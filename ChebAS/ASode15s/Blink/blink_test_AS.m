domain = [-1 1;-1 1];
cheb_struct.domain = domain;
cheb_struct.degs = [20 20];
cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

odetol = 1e-4;
tspan = [0 1];

pctClosed = 0.75;

BoundaryH = 13;

%Test with 4 patches
 Tree = ChebPatch(cheb_struct);
 Tree = Tree.split(1);
 Tree.split(2);
Tree.split(1);
 Tree.split(2);
 
H = PUchebfun(Tree);

%H = PUchebfun(@(x,y)exp(-x.^20./(1-x.^20)).*exp(-y.^20./(1-y.^20)),[-1 1;-1 1],'Degree',[20 20],'CoarseDegree',[9 9],'tol',1e-3);
%H.reset();

H.sample(@(x,y) zeros(size(x)));

P = H.copy();

setInterpMatrices(H,false);
setInterpMatrices(P,false);

[Blinks,M,y0] = SetBlinks(H,P,pctClosed,BoundaryH);

%[y,yp] = GetInitialSlope(M,y0,zeros(length(y0),1),0,{H,P},Blinks,1e-3);

 opt = odeset('mass',M,'reltol',odetol,'abstol',odetol);
 [t,U] = ASode15s(true,Blinks,tspan,y0,{H,P},opt);
