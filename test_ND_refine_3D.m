c = 0.3;
f = @(x) atan((x(:,1)+x(:,2)+x(:,3))/c);
f2 = @(x,y,z) atan((x+y+z)/c);

%f = @(x) atan((x(:,1)+x(:,2).^2)/c);
%f2 = @(x,y,z) atan((x+y.^2)/c);

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
% E = abs(ef-f2(X,Y,Z));
% % 
% max(E(:))
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