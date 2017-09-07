function [ output ] = LaplacianForward(Tree,dim,domain,sol)

Tree.sample(sol);

points1 = Tree.children{1}.leafGrids();
points2 = Tree.children{2}.leafGrids();

vx1 = Tree.evalfGrid(points1,1,2);
vy1 = Tree.evalfGrid(points1,2,2);

vx2 = Tree.evalfGrid(points2,1,2);
vy2 = Tree.evalfGrid(points2,2,2);

lap1 = vx1(:,:,3)+vy1(:,:,3);
lap2 = vx2(:,:,3)+vy2(:,:,3);

lap1 = lap1(:);
lap2 = lap2(:);

out_border1 = FindBoundaryIndex2D(dim,Tree.children{1}.domain(),domain);
lap1(out_border1) = sol(out_border1);

out_border2 = FindBoundaryIndex2D(dim,Tree.children{2}.domain(),domain);
sol2 = sol(prod(dim)+1:end);
lap2(out_border2) = sol2(out_border2);
output = [lap1;lap2];
end


