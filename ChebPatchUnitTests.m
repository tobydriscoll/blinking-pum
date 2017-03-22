clear;

tol = 1e-14;

%Test 2D

P = ChebPatch([-1 1;-1 1],[5 5]);
X = P.points();

f = @(x) x(:,1).^2+x(:,1).*x(:,2).^2;

fdx = @(x) 2*x(:,1)+x(:,2).^2;
fdy = @(x) 2*x(:,1).*x(:,2);

P.sample(f(X));

g = [linspace(-1,1,50)' linspace(-1,1,50)'];

v = P.evalf(g,1,0);
assert(norm(v(:,end)-f(g),inf)<tol)
% 
v = P.evalf(g,1,1);
assert(norm(v(:,end)-fdx(g),inf)<tol)
% 
v = P.evalf(g,2,1);
assert(norm(v(:,end)-fdy(g),inf)<tol)
% 
% % 
% % %Test 3D
P = ChebPatch([-1 1;-1 1;-1 1],[5 5 5]);
X = P.points();

f = @(x) x(:,1).^2.*x(:,3)+(x(:,3).^3).*(x(:,1)).*x(:,2).^2;

fdx = @(x) 2.*x(:,3).*x(:,1)+(x(:,3).^3).*x(:,2).^2;

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
% 
% %Test 4D
P = ChebPatch([-1 1;-1 1;-1 1;-1 1],[5 5 5 5]);
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


%Test 2D Split
P = ChebPatch([-1 1;-1 1],[65 65]);
X = P.points();
f = @(x) cos(pi*x(:,1))+sin(30*pi*x(:,2));
P.sample(f(X));

PU = P.splitleaf(0.1);
assert(isa(PU,'PUPatch'));
assert(PU.splitting_dim==2);

P = ChebPatch([-1 1;-1 1],[65 65]);
X = P.points();
f = @(x) cos(pi*x(:,1))+sin(pi*x(:,2));
P.sample(f(X));
PU = P.splitleaf(0.1);
assert(isa(PU,'ChebPatch'));


%Test 3D Split
P = ChebPatch([-1 1;-1 1;-1 1],[65 65 65]);
X = P.points();
f = @(x) cos(pi*x(:,1))+sin(30*pi*x(:,2))+cos(2*pi*x(:,3));
P.sample(f(X));

PU = P.splitleaf(0.1);
assert(isa(PU,'PUPatch'));
assert(PU.splitting_dim==2);

P = ChebPatch([-1 1;-1 1;-1 1],[65 65 65]);
X = P.points();
f = @(x) cos(pi*x(:,1))+sin(pi*x(:,2))+cos(pi*x(:,3));
P.sample(f(X));

PU = P.splitleaf(0.1);
assert(isa(PU,'ChebPatch'));



%Test 4D Split
%Use a smaller max degree.
P = ChebPatch([-1 1;-1 1;-1 1;-1 1],[33 33 33 33],5);
X = P.points();
f = @(x) cos(pi*x(:,1))+sin(30*pi*x(:,2))+cos(2*pi*x(:,3))+cos(2*pi*x(:,4));
P.sample(f(X));

PU = P.splitleaf(0.1);
assert(isa(PU,'PUPatch'));
assert(PU.splitting_dim==2);

P = ChebPatch([-1 1;-1 1;-1 1;-1 1],[33 33 33 33],5);
X = P.points();
f = @(x) cos(pi*x(:,1))+sin(pi*x(:,2))+cos(pi*x(:,3))+cos(pi*x(:,4));
P.sample(f(X));

PU = P.splitleaf(0.1);
assert(isa(PU,'ChebPatch'));