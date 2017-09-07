clear;

tol = 1e-14;
diff_tol = 1e-13;
%Test 2D

P = ChebPatch([-1 2;-1 2],[2 2],tol);
X = P.points();

f = @(x) x(:,1).^2+x(:,1).*x(:,2).^2;

fdx = @(x) 2*x(:,1)+x(:,2).^2;
fdy = @(x) 2*x(:,1).*x(:,2);

fdxx = @(x) 2*ones(size(x(:,1)));

P.sample(f(X));

%Test the point list
g = [linspace(-1,2,50)' linspace(-1,2,50)'];

v = P.evalf(g,1,0);
assert(norm(v(:,end)-f(g),inf)<tol);

v = P.evalf(g,1,1);
assert(norm(v(:,1)-f(g),inf)<tol)
assert(norm(v(:,end)-fdx(g),inf)<tol)

v = P.evalf(g,1,2);
assert(norm(v(:,1)-f(g),inf)<tol)
assert(norm(v(:,2)-fdx(g),inf)<tol)
assert(norm(v(:,3)-fdxx(g),inf)<tol)

v = P.evalf(g,2,1);
assert(norm(v(:,end)-fdy(g),inf)<tol)

g = {linspace(-1,2,7)',linspace(-1,2,7)'};


f = @(x,y) x.^2+x.*y.^2;

fdx = @(x,y) 2*x+y.^2;
fdy = @(x,y) 2*x.*y;

[X1,X2] = ndgrid(g{:});

v = P.evalfGrid(g,1,0);
assert(norm(v-f(X1,X2),inf)<tol);

v = P.evalfGrid(g,1,1);
assert(norm(v(:,:,1)-f(X1,X2),inf)<tol)
assert(norm(v(:,:,end)-fdx(X1,X2),inf)<diff_tol)

v = P.evalfGrid(g,2,1);
assert(norm(v(:,:,end)-fdy(X1,X2),inf)<diff_tol)

%Test evaluating on grid


% 
% % 
% % %Test 3D
P = ChebPatch([-1 1;-1 1;-1 1],[2 2 2]);
X = P.points();

f = @(x) x(:,1).^2.*x(:,3)+(x(:,3).^3).*(x(:,1)).*x(:,2).^2;

fdx = @(x) 2.*x(:,3).*x(:,1)+(x(:,3).^3).*x(:,2).^2;

fdxx = @(x) 2.*x(:,3);

fdy = @(x) 2*(x(:,3).^3).*(x(:,1)).*x(:,2);

fdz = @(x) x(:,1).^2+3*(x(:,3).^2).*(x(:,1)).*x(:,2).^2;

P.sample(f(X));

g = [linspace(-1,1,50)' linspace(-1,1,50)' linspace(-1,1,50)'];

v = P.evalf(g,1,0);
assert(norm(v(:,end)-f(g),inf)<tol);

v = P.evalf(g,1,1);
assert(norm(v(:,end)-fdx(g),inf)<tol);

v = P.evalf(g,2,1);
assert(norm(v(:,end)-fdy(g),inf)<tol);

v = P.evalf(g,3,1);
assert(norm(v(:,end)-fdz(g),inf)<tol);

v = P.evalf(g,1,2);
assert(norm(v(:,end)-fdxx(g),inf)<tol);

g = {linspace(-1,1,3)',linspace(-1,1,3)',linspace(-1,1,3)'};

f = @(x,y,z) x.^2.*z+(z.^3).*(x).*y.^2;

fdx = @(x,y,z) 2.*z.*x+(z.^3).*y.^2;

fdxx = @(x,y,z) 2.*z;

fdy = @(x,y,z) 2*(z.^3).*(x).*y;

fdz = @(x,y,z) x.^2+3*(z.^2).*(x).*y.^2;

[X1,X2,X3] = ndgrid(g{:});

v = P.evalfGrid(g,1,0);
E = v-f(X1,X2,X3);
assert(max(abs(E(:)))<tol);

v = P.evalfGrid(g,1,1);
E = v(:,:,:,end)-fdx(X1,X2,X3);
assert(max(abs(E(:)))<tol);

v = P.evalfGrid(g,1,2);
E = v(:,:,:,end)-fdxx(X1,X2,X3);
assert(max(abs(E(:)))<tol);

v = P.evalfGrid(g,2,1);
E = v(:,:,:,end)-fdy(X1,X2,X3);
assert(max(abs(E(:)))<tol);

v = P.evalfGrid(g,3,1);
E = v(:,:,:,end)-fdz(X1,X2,X3);
assert(max(abs(E(:)))<tol);


% 
% %Test 4D
P = ChebPatch([-1 1;-1 1;-1 1;-1 1],[2 2 2 2]);
X = P.points();

f = @(x) x(:,1).^3.*x(:,4)+x(:,2).^3.*x(:,1)+x(:,3).^3.*x(:,4).^2;

fdx = @(x) 3*x(:,1).^2.*x(:,4)+x(:,2).^3;

fdy = @(x) 3*x(:,2).^2.*x(:,1);

fdz = @(x) 3*x(:,3).^2.*x(:,4).^2;

fdw = @(x) x(:,1).^3+2*x(:,3).^3.*x(:,4);

P.sample(f(X));

g = [linspace(-1,1,50)' linspace(-1,1,50)' linspace(-1,1,50)' linspace(-1,1,50)'];

v = P.evalf(g,1,0);
assert(norm(v(:,end)-f(g),inf)<tol);

v = P.evalf(g,1,1);
assert(norm(v(:,end)-fdx(g),inf)<tol);

v = P.evalf(g,2,1);
assert(norm(v(:,end)-fdy(g),inf)<tol);

v = P.evalf(g,3,1);
assert(norm(v(:,end)-fdz(g),inf)<tol);

v = P.evalf(g,4,1);
assert(norm(v(:,end)-fdw(g),inf)<tol);


g = {linspace(-1,1,3)',linspace(-1,1,3)',linspace(-1,1,3)',linspace(-1,1,3)'};

f = @(x,y,z,w) x.^3.*w+y.^3.*x+z.^3.*w.^2;

fdx = @(x,y,z,w) 3*x.^2.*w+y.^3;

fdy = @(x,y,z,w) 3*y.^2.*x;

fdz = @(x,y,z,w) 3*z.^2.*w.^2;

fdw = @(x,y,z,w) x.^3+2*z.^3.*w;


%Test 2D Split
P = ChebPatch([-1 1;-1 1]);
X = P.points();
f = @(x) cos(pi*x(:,1))+sin(30*pi*x(:,2));
P.sample(f(X));

PU = P.splitleaf();
assert(isa(PU,'PUPatch'));
assert(PU.splitting_dim==2);

P = ChebPatch([-1 1;-1 1]);
X = P.points();
f = @(x) cos(pi*x(:,1))+sin(pi*x(:,2));
P.sample(f(X));
PU = P.splitleaf();
assert(isa(PU,'ChebPatch'));


%Test 3D Split
P = ChebPatch([-1 1;-1 1;-1 1]);
X = P.points();
f = @(x) cos(pi*x(:,1))+sin(30*pi*x(:,2))+cos(2*pi*x(:,3));
P.sample(f(X));

PU = P.splitleaf();
assert(isa(PU,'PUPatch'));
assert(PU.splitting_dim==2);

P = ChebPatch([-1 1;-1 1;-1 1]);
X = P.points();
f = @(x) cos(pi*x(:,1))+sin(pi*x(:,2))+cos(pi*x(:,3));
P.sample(f(X));

PU = P.splitleaf();
assert(isa(PU,'ChebPatch'));



%Test 4D Split
%Use a smaller max degree.
P = ChebPatch([-1 1;-1 1;-1 1;-1 1],[5 5 5 5]);
X = P.points();
f = @(x) cos(pi*x(:,1))+sin(30*pi*x(:,2))+cos(2*pi*x(:,3))+cos(2*pi*x(:,4));
P.sample(f(X));

PU = P.splitleaf();
assert(isa(PU,'PUPatch'));
assert(PU.splitting_dim==2);

% P = ChebPatch([-1 1;-1 1;-1 1;-1 1],[5 5 5 5]);
% X = P.points();
% f = @(x) cos(pi*x(:,1))+sin(pi*x(:,2))+cos(pi*x(:,3))+cos(pi*x(:,4));
% P.sample(f(X));
% 
% PU = P.splitleaf();
% assert(isa(PU,'ChebPatch'));