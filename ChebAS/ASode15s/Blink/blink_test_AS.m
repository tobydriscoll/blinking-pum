domain = [-1 1;-1 1];
cheb_struct.domain = domain;

cheb_struct.degs = [25 25];

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
flux_in_out = 4;

%Test with 4 patches
Tree = ChebPatch(cheb_struct);

overlap = 0.2;

   %This replaces Tree with one leaf overlap away from a boundary
   %| |      |
   %| |      |
   %| |      |
   %| |      |
   Tree = Tree.split(2,false,overlap);
   
   %This does the same but with the larger leaf
   %| |     | |
   %| |     | |
   %| |     | |
   %| |     | |
   Tree.children{2} = Tree.children{2}.split(2,false,overlap);
   
   
   %This does the same but with the larger leaf
   %| |     | |
   %| |     | |
   %| |_ _ _| |
   %| |     | |
   Tree.children{2}.children{1} = Tree.children{2}.children{1}.split(1,false,overlap);
   
   
   %This does the same but with the larger left over leaf
   %| |_ _ _| |
   %| |     | |
   %| |_ _ _| |
   %| |     | |
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

tspan = [0 0.733];


opt = odeset('mass',M,'reltol',odetol,'abstol',odetol);

[t,U] = ASode15s(true,Blinks,tspan,y0,{H,P},1,opt);
%save('~/Dropbox/results_small_00_no_flux.mat','Blinks','H','P','t','U');
% %  %% 

 

