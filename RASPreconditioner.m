function [ output ] = RASPreconditioner(Tree,dim,domain,cheb_length,sol)

Tree.sample(sol);
points1 = Tree.children{1}.points;
points2 = Tree.children{2}.points;

force = @(x) zeros(length(x),1);

border = @(x) exp(x(:,1)).*sin(x(:,2));

%figure out the border index
[out_border_1,in_border_1] = FindBoundaryIndex2D(dim,Tree.children{1}.domain(),domain);
[out_border_2,in_border_2] = FindBoundaryIndex2D(dim,Tree.children{2}.domain(),domain);

Dxx = kron(eye(cheb_length),diffmat(cheb_length,2));
Dyy = kron(diffmat(cheb_length,2),eye(cheb_length));

%Find the interior boundary
bval1 = Tree.evalf(points1(in_border_1,:),1,0);
bval2 = Tree.evalf(points2(in_border_2,:),1,0);

outbval1 = border(points1(out_border_1,:));
outbval2 = border(points2(out_border_2,:));

scale1x = 2/diff(Tree.children{1}.domain(1,:));
scale2x = 2/diff(Tree.children{2}.domain(1,:));

lap1 = scale1x^2*Dxx+Dyy;
E1 = eye(prod(dim));
A1 = [lap1(~(in_border_1 | out_border_1),:);E1(in_border_1,:);E1(out_border_1,:)];
b1 = [force(points1(~in_border_1 & ~out_border_1,:));bval1;outbval1];

sol1 = A1\b1;

lap2 = scale2x*Dxx+Dyy;
E2 = eye(prod(dim));
A2 = [lap2(~(in_border_2 | out_border_2),:);E2(in_border_2,:);E2(out_border_2,:)];
b2 = [force(points2(~in_border_2 & ~out_border_2,:));bval2;outbval2];

sol2 = A2\b2;

output = [sol1;sol2];
end

