NT = 1;

c = 5e-2;

%f = @(x) atan((x(:,1).^2+x(:,2)+x(:,3).^2)/c);
%f2 = @(x,y,z) atan((x.^2+y+z.^2)/c);

%f = @(x) log(1+(x(:,1).^2+x(:,2)+x(:,3).^2)/c);
%f2 = @(x,y,z) log(1+(x.^2+y+z.^2)/c);

%f = @(x) atan((x(:,1)+x(:,2).^2+x(:,3))/c);
%f2 = @(x,y,z) atan((x+y.^2+z)/c);

%f = @(x) atan((x(:,1)+x(:,2).^2)/c);
%f2 = @(x,y,z) atan((x+y.^2)/c);

%f = @(x) 1./(c+x(:,1).^2+x(:,2).^2+x(:,3).^2);
%f2 = @(x,y,z) 1./(c+x.^2+y.^2+z.^2);

%f = @(x) 1./(1+x(:,1).^2+x(:,2).^2+x(:,3).^2);
%f2 = @(x,y,z) 1./(1+x.^2+y.^2+z.^2);

u = [0.75 0.25 -0.75];

a = [15 15 15];

%f2 = @(x,y,z) cos(2*pi*u(1)+a(1)*x+a(2)*y+a(3)*z);
%f = @(x) f2(x(:,1),x(:,2),x(:,3));

f2 = @(x,y,z)  1./(a(1)^-2+(x-u(1)).^2).*1./(a(2)^-2+(y-u(2)).^2).*1./(a(3)^-2+(z-u(3)).^2);
f = @(x) f2(x(:,1),x(:,2),x(:,3));
 
%f2 = @(x,y,z) 1./(1+a(1)*x+a(2)*y+a(3)*z).^(-3);
%f = @(x) f2(x(:,1),x(:,2),x(:,3));
% 
 %f2 = @(x,y,z) exp(-a(1)*(x-u(1)).^2-a(2)*(y-u(2)).^2-a(3)*(z-u(3)).^2);
 %f = @(x) f2(x(:,1),x(:,2),x(:,3));

%f = @(x) cos(2*pi*x(:,1)).^2+cos(2*pi*x(:,2)).^2+cos(2*pi*x(:,3)).^2;
%f2 = @(x,y,z) cos(2*pi*x).^2+cos(2*pi*y).^2+cos(2*pi*z).^2;

%f2 = @(x,y,z) exp(sin(50*x))+sin(60*exp(y))+sin(70*sin(x))+sin(sin(80*y))-sin(10*(x+y))+(x.^2+y.^2)/4;
%f = @(x) f2(x(:,1),x(:,2),x(:,3));

%f2 = @(x,y,z) exp(sin(50*x))+sin(60*exp(y)).*sin(60*z)+sin(70*sin(x)).*cos(10*z)+sin(sin(80*y))-sin(10*(x+z))+(x.^2+y.^2+z.^2)/4;
%f = @(x) f2(x(:,1),x(:,2),x(:,3));


%f2 = @(x,y,z) x.*z + x.^2.*y;
%f = @(x) f2(x(:,1),x(:,2),x(:,3));

%f2 = @(x,y,z) cos(50*pi*(x+y+z));
%f = @(x) f2(x(:,1),x(:,2),x(:,3));

%f2 = @(x,y,z) 1./(1+25*(x.^2+y.^2+z.^2));
%f = @(x) f2(x(:,1),x(:,2),x(:,3));

%f2 = @(x,y,z) tanh(5*(x+z)).*exp(y);
%f = @(x) f2(x(:,1),x(:,2),x(:,3));

%f2 = @(x,y,z) 1./cosh(3*(x+y+z)).^2;
%f = @(x) f2(x(:,1),x(:,2),x(:,3));

%f2 = @(x,y,z) atan((x+y.^2+z.^2)/1e-1);
%f = @(x) f2(x(:,1),x(:,2),x(:,3));

%f2 = @(x,y,z) atan((x+y+z)*3.5);
%f = @(x) f2(x(:,1),x(:,2),x(:,3));

f2 = @(x,y,z) 1./(1+25*(x.^2+y.^2+z.^2));
f = @(x) f2(x(:,1),x(:,2),x(:,3));

TIMES = zeros(NT,1);

for i=1:NT
tic;
TREE = PUFun([-1 1;-1 1;-1 1],[6 6 6],f,1e-12);
TIMES(i) = toc;
end
% 
mean(TIMES)

x = linspace(-1,1,50)';

G = {x x x};

TIMES = zeros(NT,1);

for i=1:NT
tic;
ef = TREE.ChebRoot.evalfGrid(G);
TIMES(i) = toc;
end

mean(TIMES)
% % 
[X,Y,Z] = ndgrid(x,x,x);
% 
% 
% surf(X,Y,ef);
% 
% figure();
% TREE.ChebRoot.plotdomain();
E = abs(ef-f2(X,Y,Z));
% 
max(E(:))

length(TREE)

% for i=1:NT
% tic;
%  F = chebfun3(f2);
%  TIMES(i) = toc;
% end
% 
% mean(TIMES)
% % % % 
% for i=1:NT
% tic;
% ef2 = F(X,Y,Z);
% TIMES(i) =toc;
% end

%mean(TIMES)

%E2 = abs(ef2 - f2(X,Y,Z));

%max(E2(:))

%rank(F)