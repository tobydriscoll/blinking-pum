clear;

x=chebpts(3);
y=chebpts(5);
z=chebpts(17);
w=chebpts(33);

for i=1:3
    C{i}=x;
end

xc = linspace(-1,1,6)';
yc = linspace(-1,1,6)';
zc = linspace(-1,1,6)';
wc = linspace(-1,1,6)';

XP = [xc yc];

Mx = barymat(xc,x);
My = barymat(yc,y);
Mz = barymat(zc,z);

[X,Y] = ndgrid(x,y);

F = X.^2+Y.*X;

for j=1:6
G = F;
%G = reshape(feval(chebfun(G),XP(j,1)),[length(y) length(z) length(w)]);
%G = reshape(feval(diff(chebfun(G),1),XP(j,1)),[length(y) length(z)]);
G = feval(diff(chebfun(G),1),XP(j,1));
EV(j) = feval(diff(chebfun(G'),0),XP(j,2));
end

EV'

2*xc+yc



