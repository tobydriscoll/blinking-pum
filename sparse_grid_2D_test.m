domain = [-1 1;-1 1];

u = [0.75 0.25];

a = [5 10];

tol = 1e-12;

opt = spset('GridType','Chebyshev','reltol',tol,'abstol',tol,...
    'maxdepth',10,'vectorized','on');
    
test_funs = {@(x,y) log(1+10^5*(x.^2+y.^2));
            @(x,y) atan(10^2*(x+y.^2));
            @(x,y) 10^(-4)./((10^(-4)+x.^2).*(10^(-4)+y.^2));
            @(x,y) franke(x,y);
            @(x,y) cos(2*pi*u(1)+a(1)*x+a(2)*y);
            @(x,y)  1./(a(1)^-2+(x-u(1)).^2).*1./(a(2)^-2+(y-u(2)).^2);
            @(x,y) 1./(1+a(1)*x+a(2)*y).^(-3);
            @(x,y) exp(-a(1)*(x-u(1)).^2-a(2)*(y-u(2)).^2)};
        
CON_TIME = zeros(length(test_funs),1);
INTERP_TIME = zeros(length(test_funs),1);

for i=1:length(test_funs)
    tic,SP = spvals(test_funs{i},2,domain,opt); INTERP_TIME(i) = toc;    
    
end