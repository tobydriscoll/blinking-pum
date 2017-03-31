%For simplicity, the function should be written like this
%I will do something clever later.
b = 0.05;
f = @(x) atan((x(:,1)+x(:,2))/b);
func = @(x,y) atan((x+y)/b);

degs_in = [6 6];
domain = [-1 1;-1 1];
tic
PUf = PUFun(domain,degs_in,f);
toc

tic
fc=chebfun2(func)
toc


f = @(x) atan((x(:,1)/b)+atan(x(:,2))/b);
func = @(x,y) atan(x/b)+atan(y/b);

tic
PUf = PUFun(domain,degs_in,f);
toc

tic
fc=chebfun2(func)
toc