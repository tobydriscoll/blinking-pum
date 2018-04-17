%DOMAIN = Disk(1,[0 0]);
DOMAIN = DoubleAstroid();
%DOMAIN = Diamond();

delta = 0;
OUTERBOX = [-0.95 0.95;-0.95 0.95];
%f = @(x,y)1./(((x-1.1).^2)+(y-1.1).^2).^2;
f =  @(x,y) atan(3*(x.^2+y));
%f2 = @(x,y)1./(((x+1.05).^2)+(y+1.05).^2).^2;
%f = @(x,y) exp(x+y);
%f = @(x,y) cos(24*x - 32*y).*sin(21*x - 28*y);
%f = @(x,y) abs(x.*y);
% fx = @(x,y) (2*x)./(-1.1 + x.^2 + y.^2);
% 
% DOMAIN = SquareSlice(-0.5,[-1 1;-1 1]);
% OUTERBOX = [-1 1;-1 1];
% f = @(x,y) atan((x+y-0.5)*5);

% DOMAIN = Heart();
% OUTERBOX = [-3.2,3.2;0,4.7];
% f = @(x,y) atan((x+y-3.1)*10);

%DOMAIN = Ball(1,[0 0 0]);
%OUTERBOX = [-1 1;-1 1;-1 1];
%f = @(x,y,z) atan(3*(x+y+z));

tic,TREELS = PUFunLS(f,DOMAIN,OUTERBOX,'degreeIndex',[4 4],'ChebDegreeIndex',[6 6],'tol',1e-10);toc
%tic,TREELS1 = PUFunLS(f1,DOMAIN,OUTERBOX,'degreeIndex',[4 4],'ChebDegreeIndex',[6 6],'tol',1e-12);toc
%tic,TREELS2 = PUFunLS(f2,DOMAIN,OUTERBOX,'degreeIndex',[4 4],'ChebDegreeIndex',[6 6],'tol',1e-12);toc


%tic, TREE3 = TREELS1/TREELS2; toc
%plot(DOMAIN);hold on; plotdomain(TREELS.ChebRoot); hold off;

%x = OUTERBOX(1,1) + diff(OUTERBOX(1,:))*rand(200,1);
%y = OUTERBOX(2,1) + diff(OUTERBOX(2,:))*rand(200,1);

x = linspace(-1,1,100)';
y = linspace(-1,1,100)';


% 
[X,Y] = ndgrid(x,y);
% 
V = TREELS.ChebRoot.evalfGrid({x y});

F = f(X,Y);

E = V - F;

ind = DOMAIN.Interior([X(:) Y(:)]);

B = DOMAIN.Boundary(800);

VB = TREELS.ChebRoot.evalf(B);

NP = sum(ind);

P = [X(ind) Y(ind)];

TRI = delaunayTriangulation([B;P],[(1:length(B)-1)' (2:length(B))'; length(B) 1]);

TF = isInterior(TRI);

max(abs(E(ind)))/max(abs(F(ind)))

trisurf(TRI.ConnectivityList(TF,:),TRI.Points(:,1),TRI.Points(:,2),[VB;V(ind)]);

% DiffTree = TREELS.diff(1,1);
% 
% V2 = DiffTree.ChebRoot.evalfGrid({x y});
% 
% V2 = V2 - fx(X,Y);
% 
% max(abs(V2(ind)))