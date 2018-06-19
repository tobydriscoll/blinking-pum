DOMAIN = {Disk(1,[0 0]) , Diamond(), DoubleAstroid()};
OUTERBOX = {[-1 1;-1 1], [-0.95 0.95;-0.95 0.95] , [-1 1;-1 1]};
funs = {@(x,y) exp(x+y), @(x,y)1./(((x-1.1).^2)+(y-1.1).^2).^2, @(x,y) cos(24*x-32*y).*sin(21*x-28*y),@(x,y) atan(3*(x.^2+y))};

CONS_TIMES = zeros(length(DOMAIN),length(funs));
INTERP_TIMES = zeros(length(DOMAIN), length(funs));
INTERP_ERR = zeros(length(DOMAIN),length(funs));
NUMPTS = zeros(length(DOMAIN),length(funs));

DegreeIND = [4 4];
ChebIND = [6 6];
tol = 1e-10;

for D_i=1:length(DOMAIN)
        for f_i=1:length(funs)
            
             tic,F = PUFunLS(funs{f_i},DOMAIN{D_i},OUTERBOX{D_i},'degreeIndex',DegreeIND,'ChebDegreeIndex',ChebIND,'tol',tol);CONS_TIMES(D_i,f_i) = toc;
             
             x = linspace(OUTERBOX{D_i}(1,1),OUTERBOX{D_i}(1,2),200)';
             y = linspace(OUTERBOX{D_i}(2,1),OUTERBOX{D_i}(2,2),200)';
             [X,Y] = ndgrid(x,y);
             
             tic, V = F.evalfGrid({x,y}); INTERP_TIMES(D_i,f_i) = toc;
             IN_IND = DOMAIN{D_i}.Interior([X(:) Y(:)]);
             
             E = V - funs{f_i}(X,Y);
             
             INTERP_ERR(D_i,f_i) = max(abs(E(IN_IND)));
             
             NUMPTS(D_i,f_i) = length(F);
        end
end

STATS = [INTERP_ERR(:,1) CONS_TIMES(:,1) INTERP_TIMES(:,1) NUMPTS(:,1);
         INTERP_ERR(:,2) CONS_TIMES(:,2) INTERP_TIMES(:,2) NUMPTS(:,2);
         INTERP_ERR(:,3) CONS_TIMES(:,3) INTERP_TIMES(:,3) NUMPTS(:,3);
         INTERP_ERR(:,4) CONS_TIMES(:,4) INTERP_TIMES(:,4) NUMPTS(:,4)];


csvwrite('LS2D_stats_2.csv',STATS);