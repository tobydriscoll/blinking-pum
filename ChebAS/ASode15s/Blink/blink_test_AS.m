domain = [-1 1;-1 1];
cheb_struct.domain = domain;

cheb_struct.degs = [20 20];

cheb_struct.cdegs = [9 9];
cheb_struct.split_flag = [true true];
cheb_struct.tol = 1e-4;

odetol = 1e-4;
tspan = [0 0.3];

pctClosed = 0.7;

pA = 6.11e-6;
pS = 3.09e-6;

he = 2;


BoundaryH = 13;
initial_volume = 24;
flux_in_out = 1;

%Test with 4 patches
Tree = ChebPatch(cheb_struct);

overlap = 0.2;

% % 
%Tree = Tree.split(1);
%Tree.split(2);
%Tree.split(1);
%Tree.split(2);
% Tree.split(1);
% Tree.split(2);
% 
   Tree = Tree.split(2,false,overlap);
   Tree.children{2} = Tree.children{2}.split(2,false,overlap);
   Tree.children{2}.children{1} = Tree.children{2}.children{1}.split(1,false,overlap);
   Tree.children{2}.children{1}.children{2} = Tree.children{2}.children{1}.children{2}.split(1,false,overlap);

  Tree.clean();

H = PUchebfun(Tree);

%H = PUchebfun(@(x,y)exp(-x.^20./(1-x.^20)).*exp(-y.^20./(1-y.^20)),[-1 1;-1 1],'Degree',[20 20],'CoarseDegree',[9 9],'tol',1e-3);
%H.reset();

H.sample(@(x,y) zeros(size(x)));

P = H.copy();

setInterpMatrices(H,false);
setInterpMatrices(P,false);

[Blinks,M,y0] = setBlinks(H,P,pctClosed,BoundaryH,pA,pS,he,initial_volume,flux_in_out);
% 
% %[y,yp] = GetInitialSlope(M,y0,zeros(length(y0),1),0,{H,P},Blinks,1e-3);
% 
tspan = [0 Blinks{1}.period];
%tspan = [0 1e-9];

opt = odeset('mass',M,'reltol',odetol,'abstol',odetol);

[t,U] = ASode15s(false,Blinks,tspan,y0,{H,P},1,opt);
save('~/Dropbox/results_small_07.mat','Blinks','H','P','t','U');
% %  %% 

 

