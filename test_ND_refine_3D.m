c = 0.1;
f = @(x) atan((x(:,1)+x(:,2)+x(:,3))/c);
f2 = @(x,y,z) atan((x+y+z)/c);

%f = @(x) atan((x(:,1)+x(:,2).^2)/c);
%f2 = @(x,y,z) atan((x+y.^2)/c);

%f = @(x) 1./(1+x(:,1).^2+x(:,2).^2+x(:,3).^2);
%f2 = @(x,y,z) 1./(1+x.^2+y.^2+z.^2);

tic;
TREE = PUFun([-1 1;-1 1;-1 1],[5 5 5],f,1e-10);
toc;

x = linspace(-1,1,100)';
% 
G = {x x x};
% 
tic;
ef = TREE.evalfGrid(G,1,0);
toc
% 
[X,Y,Z] = ndgrid(x,x,x);

%surf(X,Y,ef);

%figure();
%TREE.ChebRoot.plotdomain();
% % 
E = abs(ef-f2(X,Y,Z));
% % 
max(E(:))
% tic;
% F = chebfun3t(f2,'eps',1e-10)
% toc;
% % 
% tic;
% ef2 = F(X,Y);
% toc;
% 
% E2 = abs(ef2-f2(X,Y));
% 
% max(E2(:))