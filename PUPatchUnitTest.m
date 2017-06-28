clear;

tol = 1e-14;

diff_tol = 1e-11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 2D Split
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = ChebPatch([-1 1;-1 1],[6 6],[true true],tol);
X = P.points();
b = 0.5;
f = @(x) x(:,1).*atan(x(:,2)./b);
fx = @(x) atan(x(:,2)./b);
fy = @(x) b*x(:,1)./(b^2+x(:,2).^2);
P.sample(f(X));

PU = P.splitleaf();
%This should split
assert(isa(PU,'PUPatch'));
assert(PU.splitting_dim == 2);

X = PU.points();
PU.sample(f(X));

PU.PUsplit();

%The leaves after the split should be refined
assert(PU.is_refined);
assert(isa(PU.children{1},'ChebPatch'));
assert(isa(PU.children{2},'ChebPatch'));
% 
X = PU.points();
% 
v = PU.evalf(X,1,1);
assert(norm(v(:,1)-f(X),inf)<tol);
assert(norm(v(:,2)-fx(X),inf)<diff_tol);

v = PU.evalf(X,2,1);
assert(norm(v(:,2)-fy(X),inf)<diff_tol);

g = {linspace(-1,1,60)',linspace(-1,1,60)'};
[X1,X2] = ndgrid(g{:});

f = @(x,y) x.*atan(y./b);
fx = @(x,y) atan(y./b);
fy = @(x,y) b*x./(b^2+y.^2);

v = PU.evalfGrid(g,1,1);

E=v(:,:,1)-f(X1,X2);
assert(max(abs(E(:)))<tol);

E=v(:,:,2)-fx(X1,X2);
assert(max(abs(E(:)))<diff_tol);

v = PU.evalfGrid(g,2,1);
E=v(:,:,2)-fy(X1,X2);
assert(max(abs(E(:)))<diff_tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 3D Split
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = ChebPatch([-1 1;-1 1;-1 1],[6 6 6],[true true true],tol);
X = P.points();
b = 0.5;
f = @(x) x(:,1).*atan(x(:,2)./b)+x(:,3).^6;
fx = @(x) atan(x(:,2)./b);
fy = @(x) b*x(:,1)./(b^2+x(:,2).^2);
fz = @(x) 6*x(:,3).^5;
P.sample(f(X));

PU = P.splitleaf();
%This should split
assert(isa(PU,'PUPatch'));
assert(PU.splitting_dim == 2);

X = PU.points();
PU.sample(f(X));

PU.PUsplit();

%The leaves after the split should be refined
assert(PU.is_refined);
assert(isa(PU.children{1},'ChebPatch'));
assert(isa(PU.children{2},'ChebPatch'));
% 
X = PU.points();
% 
v = PU.evalf(X,1,1);
assert(norm(v(:,1)-f(X),inf)<tol);
assert(norm(v(:,2)-fx(X),inf)<diff_tol);
% 
v = PU.evalf(X,2,1);
assert(norm(v(:,2)-fy(X),inf)<diff_tol);

v = PU.evalf(X,3,1);

assert(norm(v(:,2)-fz(X),inf)<diff_tol);

g = {linspace(-1,1,60)',linspace(-1,1,60)',linspace(-1,1,60)'};

f = @(x,y,z) x.*atan(y./b)+z.^6;
fx = @(x,y,z) atan(y./b);
fy = @(x,y,z) b*x./(b^2+y.^2);
fz = @(x,y,z) 6*z.^5;

[X1,X2,X3] = ndgrid(g{:});

v = PU.evalfGrid(g,1,1);

E=v(:,:,:,1)-f(X1,X2,X3);
assert(max(abs(E(:)))<tol);

E=v(:,:,:,2)-fx(X1,X2,X3);
assert(max(abs(E(:)))<diff_tol);

v = PU.evalfGrid(g,2,1);
E=v(:,:,:,2)-fy(X1,X2,X3);
assert(max(abs(E(:)))<diff_tol);

v = PU.evalfGrid(g,3,1);
E=v(:,:,:,2)-fz(X1,X2,X3);
assert(max(abs(E(:)))<diff_tol);



