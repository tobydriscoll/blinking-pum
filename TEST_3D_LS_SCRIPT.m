DOMAIN = {Ball(1,[0 0 0])};
OUTERBOX = {[-1 1;-1 1;-1 1]};
funs = {@(x,y,z) exp(x+y+z), @(x,y,z)atan(3*(x+y+z)), @(x,y,z) log(1.2-(x.^2+y.^2+z.^2)),@(x,y,z) 1./((x-1.1).^2+(y-1.1).^2+(z-1.1).^2).^2};

CONS_TIMES = zeros(length(DOMAIN),length(funs));
INTERP_TIMES = zeros(length(DOMAIN), length(funs));
INTERP_ERR = zeros(length(DOMAIN),length(funs));
NUMPTS = zeros(length(DOMAIN),length(funs));

DegreeIND = [3 3 3];
ChebIND = [5 5 5];
tol = 1e-3;

for D_i=1:length(DOMAIN)
        for f_i=1:length(funs)
            
             tic,F = PUFunLS(funs{f_i},DOMAIN{D_i},OUTERBOX{D_i},'degreeIndex',DegreeIND,'ChebDegreeIndex',ChebIND,'tol',tol);CONS_TIMES(D_i,f_i) = toc;
             
             x = linspace(OUTERBOX{D_i}(1,1),OUTERBOX{D_i}(1,2),200)';
             y = linspace(OUTERBOX{D_i}(2,1),OUTERBOX{D_i}(2,2),200)';
             z = linspace(OUTERBOX{D_i}(3,1),OUTERBOX{D_i}(3,2),200)';
             
             [X,Y,Z] = ndgrid(x,y,z);
             
             tic, V = F.evalfGrid({x,y,z}); INTERP_TIMES(D_i,f_i) = toc;
             IN_IND = DOMAIN{D_i}.Interior([X(:) Y(:) Z(:)]);
             
             E = V - funs{f_i}(X,Y,Z);
             
             INTERP_ERR(D_i,f_i) = max(abs(E(IN_IND)));
             
             NUMPTS(D_i,f_i) = length(F);
        end
end

STATS = [INTERP_ERR(:,1) CONS_TIMES(:,1) INTERP_TIMES(:,1) NUMPTS(:,1);
         INTERP_ERR(:,2) CONS_TIMES(:,2) INTERP_TIMES(:,2) NUMPTS(:,2);
         INTERP_ERR(:,3) CONS_TIMES(:,3) INTERP_TIMES(:,3) NUMPTS(:,3);
         INTERP_ERR(:,4) CONS_TIMES(:,4) INTERP_TIMES(:,4) NUMPTS(:,4)];


csvwrite('LS3D_stats_1.csv',STATS);