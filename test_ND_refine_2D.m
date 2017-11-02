c = 0.01;

NT = 1;

NT2 = 0;

%f = @(x) log((x(:,1).^2+x(:,2).^4)/c+1);
%f2 = @(x,y) log((x.^2+y.^4)/c+1);

%f = @(x) atan((x(:,1)+x(:,2).^2)/c);
%f2 = @(x,y) atan((x+y.^2)/c);

f = @(x) atan((x(:,1).^2+x(:,2).^2)/c);
f2 = @(x,y) atan((x.^2+y.^2)/c);

%f = @(x) exp(-100*(x(:,1).^2+x(:,2).^2));
%f2 = @(x,y) exp(-100*(x.^2+y.^2));
%f = @(x) exp(x(:,1)).*sin(x(:,2));
%f2 = @(x,y) exp(x).*sin(y);

%f = @(x) atan((x(:,1))/c);
%f2 = @(x,y) atan((x)/c);

%f = @(x) (c./(c+x(:,1).^2)).*(c./(c+x(:,2).^2));
%f2 = @(x,y) (c./(c+x.^2)).*(c./(c+y.^2));

%f = @(x) (c^2./(c^2+(x(:,1)+x(:,2)).^2)).*(c^2./(c^2+(x(:,1)-x(:,2)).^2));
%f2 = @(x,y) (c^2./(c^2+(x+y).^2)).*(c^2./(c^2+(x-y).^2));

TIMES = zeros(NT,1);

for i=1:NT
tic;
TREE = PUFun([-1 1;-1 1],[6 6],f,1e-10);
TIMES(i)=toc;
end

mean(TIMES)

x = linspace(-1,1,200)';

G = {x x};


TIMESEV = zeros(NT,1);

for i=1:NT
tic;
ef = TREE.ChebRoot.evalfZoneGrid({x x});
TIMESEV(i) = toc;
end

mean(TIMESEV)

[X,Y] = ndgrid(x,x);

defaultOpts = {'facecolor', 'flat', 'edgealpha', .5, 'edgecolor', 'none'};


surf(X,Y,ef,defaultOpts{:});

E = abs(ef-f2(X,Y));

max(E(:))

TIMES = zeros(NT2,1);

for i=1:NT2
tic;
F = chebfun2(f2);
TIMES(i) = toc;
end

mean(TIMES)

TIMESEV = zeros(NT2,1);


for i=1:NT2
tic;
ef2 = F(X,Y);
TIMESEV(i) = toc;
end

mean(TIMESEV)

E2 = abs(ef2-f2(X,Y));

max(E2(:))