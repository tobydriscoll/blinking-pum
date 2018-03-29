%DOMAIN = Disk(1,[0 0]);
%DOMAIN = Astroid(pi/4);
DOMAIN = DoubleAstroid();
%DOMAIN = ParabRegion(-0.9);
delta = 0;
OUTERBOX = [-(1+delta) 1+delta;-(1+delta) 1+delta];
f = @(x,y)1./(((x-1.1).^2)+(y-1.1).^2).^2;
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

tic,TREELS = PUFunLS(f,DOMAIN,OUTERBOX,'degreeIndex',[4 4],'ChebDegreeIndex',[6 6],'tol',1e-12);toc
plot(DOMAIN);hold on; plotdomain(TREELS.ChebRoot); hold off;

%x = OUTERBOX(1,1) + diff(OUTERBOX(1,:))*rand(200,1);
%y = OUTERBOX(2,1) + diff(OUTERBOX(2,:))*rand(200,1);

x = linspace(-1,1,200)';
y = linspace(-1,1,200)';

[X,Y] = ndgrid(x,y);
V = TREELS.ChebRoot.evalfGrid({x y});

F = f(X,Y);

V = V - f(X,Y);

ind = DOMAIN.Interior([X(:) Y(:)]);

max(abs(V(ind)))/max(abs(F(ind)))

% DiffTree = TREELS.diff(1,1);
% 
% V2 = DiffTree.ChebRoot.evalfGrid({x y});
% 
% V2 = V2 - fx(X,Y);
% 
% max(abs(V2(ind)))