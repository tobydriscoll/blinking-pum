c = 5e-2;

%f = @(x) atan((x(:,1).^2+x(:,2)+x(:,3).^2)/c);
%f2 = @(x,y,z) atan((x.^2+y+z.^2)/c);

%f = @(x) log(1+(x(:,1).^2+x(:,2)+x(:,3).^2)/c);
%f2 = @(x,y,z) log(1+(x.^2+y+z.^2)/c);

f = @(x) atan((x(:,1)+x(:,2).^2+x(:,3))/c);
f2 = @(x,y,z) atan((x+y.^2+z)/c);

%f = @(x) atan((x(:,1)+x(:,2).^2)/c);
%f2 = @(x,y,z) atan((x+y.^2)/c);

%f = @(x) 1./(c+x(:,1).^2+x(:,2).^2+x(:,3).^2);
%f2 = @(x,y,z) 1./(c+x.^2+y.^2+z.^2);

%f = @(x) 1./(1+x(:,1).^2+x(:,2).^2+x(:,3).^2);
%f2 = @(x,y,z) 1./(1+x.^2+y.^2+z.^2);

tic;
TREE = PUFun([-1 1;-1 1;-1 1],[6 6 6],f,1e-6);
toc;

x = linspace(-1,1,200)';

G = {x x x};

tic;
ef = TREE.evalfGrid(G);
toc
% % 
[X,Y,Z] = ndgrid(x,x,x);
% 
% 
% %surf(X,Y,ef);
% 
% %figure();
% %TREE.ChebRoot.plotdomain();
E = abs(ef-f2(X,Y,Z));
% 
max(E(:))

tic;
F = chebfun3(f2,'eps',1e-10);
toc;
% % 
tic;
ef2 = F(X,Y,Z);
toc;

E2 = abs(ef2 - f2(X,Y,Z));

max(E2(:))
