c = 7.5e-2;

NT = 1;

NT2 = 1;

ep = 1e-5;

degs = [7 7];

domain = [-1 1;-1 1];
%f = @(x) log((x(:,1).^2+x(:,2).^4)/c+1);
%f2 = @(x,y) log((x.^2+y.^4)/c+1);

%f = @(x) atan((x(:,1)+x(:,2).^2)/c);
%f2 = @(x,y) atan((x+y.^2)/c);

%alpha = 1/(10*pi);
%f  = @(x) sin(1./(alpha+sqrt(x(:,1).^2+x(:,2).^2)));
%f2 = @(x) sin(1./(alpha+sqrt(x.^2+y.^2)));

%f = @(x) (1-exp((x(:,1)-1)/ep)).*(1-exp((x(:,2)-1)/ep)).*cos(pi*(x(:,1)+x(:,2)));
%f2 = @(x,y) (1-exp((x-1)/ep)).*(1-exp((y-1)/ep)).*cos(pi*(x+y));

%f = @(x) atan((x(:,1)+x(:,2).^2)/c);
%f2 = @(x,y) atan((x+y.^2)/c);

f  = @(x) x(:,1).^3.*exp(3.5*pi*x(:,2));
f2 = @(x,y) x.^3.*exp(3.5*pi*y);

%f = @(x) atan(x(:,2)./x(:,1));
%f2 = @(x,y) atan(y./x);

%f2 = @(x,y) atan(alpha.*((-1).*r0+((x+(-1).*xc).^2+(y+(-1).*yc).^2).^(1/2)));
%f = @(x) f2(x(:,1),x(:,2));

%f = @(x) exp(-100*(x(:,1).^2+x(:,2).^2));
%f2 = @(x,y) exp(-100*(x.^2+y.^2));
%f = @(x) exp(x(:,1)).*sin(x(:,2));
%f2 = @(x,y) exp(x).*sin(y);

%f = @(x) atan((x(:,1))/c);
%f2 = @(x,y) atan((x)/c);

%f = @(x)  c./((c+x(:,1).^2).*(c+x(:,2).^2));
%f2 = @(x,y) c./((c+x.^2).*(c+y.^2));

%f = @(x) (c^2./(c^2+(x(:,1)+x(:,2)).^2)).*(c^2./(c^2+(x(:,1)-x(:,2)).^2));
%f2 = @(x,y) (c^2./(c^2+(x+y).^2)).*(c^2./(c^2+(x-y).^2));

%f = @(x) sqrt(2-x(:,1).^2-x(:,2).^2);
%f2 = @(x,y) sqrt(2-x.^2-y.^2);

%f = @(x) franke(x(:,1),x(:,2)); 
%f2 = @(x,y) franke(x,y);

%f = @(x) sin(pi*x(:,1)+50*pi*x(:,2));
%f2 = @(x,y) sin(pi*x+50*y)

%f = @(x) 1./(1+100*(x(:,1).^2+x(:,2).^2).^2);
%f2 = @(x,y) 1./(1+100*(x.^2+y.^2).^2);


TIMES = zeros(NT,1);

for i=1:NT
    
tic;
TREE = PUFun(domain,degs,f,1e-6);
TIMES(i)=toc;

end

mean(TIMES)

x = linspace(domain(1,1),domain(1,2),100)';
y = linspace(domain(2,1),domain(2,2),100)';

G = {x y};


TIMESEV = zeros(NT,1);

for i=1:NT
tic;
ef = TREE.evalfGrid({x x});
TIMESEV(i) = toc;
end

mean(TIMESEV)

[X,Y] = ndgrid(x,x);

defaultOpts = {'facecolor', 'flat', 'edgealpha', .5, 'edgecolor', 'none'};


surf(X,Y,ef,defaultOpts{:});

E = abs(ef-f2(X,Y));

max(E(:)./max(ef(:)))

TIMES = zeros(NT2,1);

for i=1:NT2
tic;
F = chebfun2(f2,[domain(1,:) domain(2,:)]);
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