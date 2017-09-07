overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
split_flag = [1 1];
tol = 1e-12;
cheb_length  = 33;
dim = [33 33];
domain1 = [-1 overlap;-1 1];
domain2 = [-overlap 1;-1 1];

overlap_in = [-overlap overlap];

x = linspace(-1,1,100)';
[X,Y] = ndgrid(x);

children{1} = ChebPatch(domain1,deg_in,split_flag,tol);
children{2} = ChebPatch(domain2,deg_in,split_flag,tol);

points1 = children{1}.points;
points2 = children{2}.points;

force = @(x) zeros(length(x),1);

border = @(x) exp(x(:,1)).*sin(x(:,2));

TRUE = exp(X).*sin(Y);

%figure out the border index
[out_border_1,in_border_1] = FindBoundaryIndex2D(dim,domain1,domain,points1);
[out_border_2,in_border_2] = FindBoundaryIndex2D(dim,domain2,domain,points2);


Tree = PUPatch(domain,overlap_in,2*cheb_length,children,1,[]);

Dxx = kron(eye(cheb_length),diffmat(cheb_length,2));
Dyy = kron(diffmat(cheb_length,2),eye(cheb_length));

sol1 = zeros(length(points1),1);
sol2 = zeros(length(points2),1);

Tree.sample([sol1;sol2]);
%We will solve lap u = 1 with zero boundary conditions

while true

F = Tree.evalfGrid({x x},1,0);


pause(0.5);

surf(X,Y,F);
E = abs(F-TRUE);
E = max(E(:));
title(sprintf('The error is: %g',E));
%Step2: solve the problems locally

%Find the interior boundary
bval1 = Tree.evalf(points1(in_border_1,:),1,0);
bval2 = Tree.evalf(points2(in_border_2,:),1,0);

outbval1 = border(points1(out_border_1,:));
outbval2 = border(points2(out_border_2,:));

scale1x = 2/diff(domain1(1,:));
scale2x = 2/diff(domain2(1,:));

lap1 = scale1x^2*Dxx+Dyy;
E1 = eye(prod(dim));
A1 = [lap1(~(in_border_1 | out_border_1),:);E1(in_border_1,:);E1(out_border_1,:)];
b1 = [force(points1(~in_border_1 & ~out_border_1,:));bval1;outbval1];
% 
sol1 = A1\b1;
% 
lap2 = scale2x^2*Dxx+Dyy;
E2 = eye(prod(dim));
A2 = [lap2(~(in_border_2 | out_border_2),:);E2(in_border_2,:);E2(out_border_2,:)];
b2 = [force(points2(~in_border_2 & ~out_border_2,:));bval2;outbval2];

sol2 = A2\b2;

Tree.sample([sol1;sol2]);
end
