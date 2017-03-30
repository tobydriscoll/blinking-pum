b = 1e-1;
c = 1e-2;
%For simplicity, the function should be written like this
%I will do something clever later.
f = @(x) x(:,1).*atan(x(:,2)./b)+atan(x(:,3)./c);
degs_in = [6 6 6];
domain = [-1 1;-1 1;-1 1];
PUf = PUFun(domain,degs_in,f);