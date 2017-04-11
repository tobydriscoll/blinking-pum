domain = [0,1];

x = chebpts(5);

g = cos(x.^2);

G = chebfun(g,[0,10]);