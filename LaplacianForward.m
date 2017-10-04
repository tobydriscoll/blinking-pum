function [ output ] = LaplacianForward(Tree,domain,sol)

Tree.sample(sol);

LEAVES = Tree.collectLeaves({});

output = [];

step_n = 0;

for k=1:length(LEAVES)
    
    points = LEAVES{k}.leafGrids();
    
    pointsl = LEAVES{k}.points();
    
    dim = LEAVES{k}.degs;
    
    sol_k = sol(step_n+(1:prod(dim)));
        
    [out_border, in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
    
    %vx = LEAVES{k}.evalfGrid(points,1,2);
    %vy = LEAVES{k}.evalfGrid(points,2,2);
    
    %lap = vx(:,:,3)+vy(:,:,3);
    
    lap = LEAVES{k}.linOp*sol_k;
    
    %approx = Tree.evalf(pointsl(in_border,:),1,0);
    approx = Tree.evalfZone(pointsl(in_border,:));
    
    lap(in_border) = lap(in_border)-approx;
    
    output = [output;lap];
    
    step_n = step_n + prod(dim);
    
end

end


