domain = [-1 1;-1 1];
in_domain = Disk(2,[0 0]);
degs = [33 33];
fine_degs = [65 65];

B_n = 65;

tau1 = 1e-7;

Dx = ChebDiff(degs(1));
Dy = ChebDiff(degs(2));

DX = kron(eye(degs(2)),Dx);
DXY = kron(Dy,Dx);
DXX = kron(eye(degs(2)),Dx^2);
DYY = kron(Dy^2,eye(degs(1)));

tol = 1e-14;

x = chebpts(degs(1));
y = chebpts(degs(2));
[X,Y] = ndgrid(x,y);

x_f = chebpts(fine_degs(1))';
y_f = chebpts(fine_degs(2))';

[Xf,Yf] = ndgrid(x_f,y_f);

XPf = [Xf(:) Yf(:)];

grid_sq_ind = (XPf(:,1)-1).^2+(XPf(:,2)-1).^2<=4;

%fine grid in domain
XPf = XPf(grid_sq_ind,:);

theta = linspace(pi,3*pi/2,B_n)';

% boundry of disk with radius 2 centered at (1,1) 
% intersected with square.
B = [2*sin(theta)+1 2*cos(theta)+1]; 
XPf = [XPf;B];

Mx = clenshaw(XPf(:,1),eye(degs(1)));
My = clenshaw(XPf(:,2),eye(degs(2)));

%Construct Interp matrix for boundary points
M = zeros(length(XPf),prod(degs));

for i=1:length(XPf)
M(i,:) = kron(My(i,:),Mx(i,:));
end

grid_sq_ind_b =  XPf(:,1)==1 | XPf(:,2)==1 |  ((XPf(:,1)-1).^2+(XPf(:,2)-1).^2==4);

grid_sq_ind_in = ~grid_sq_ind_b;
num_in_circ = sum(grid_sq_ind);

%This is a simply matrix that
% 1) restricts the points in the grid inside
%    the domain to there function values
% 2) enforces values at the boundary
A = zeros(length(XPf),prod(degs));

A(grid_sq_ind_in,:) = M(grid_sq_ind_in,:)*(DXX+DYY);
A(grid_sq_ind_b,:) = M(grid_sq_ind_b,:);


f = @(x,y) zeros(size(x));
bound = @(x,y) log((x+1).^2+(y+1).^2);
fdx = @(x,y) 2*(1+x)./((1+x).^2+(1+y).^2);

b = zeros(length(XPf),1);
b(grid_sq_ind_in,:) = f(XPf(grid_sq_ind_in,1),XPf(grid_sq_ind_in,2));
b(grid_sq_ind_b,:) = bound(XPf(grid_sq_ind_b,1),XPf(grid_sq_ind_b,2));

bf = bound(XPf(:,1),XPf(:,2));
bfdx = fdx(XPf(:,1),XPf(:,2));

lap = DXX+DYY+2*DXY;

%tikhonov regularization
A2 = [tau1*lap;A];
b2 = [zeros(prod(degs),1);b];
tic;
V = A2\b2;
toc

%basic solution
% tic;
% V = A\b;
% toc

norm(A*V-b,inf)
norm(M*V-bf,inf)
norm(M*DX*V-bfdx,inf)

V1 = chebfun2.coeffs2vals(reshape(V,degs));

surf(X,Y,V1);
hold on;
scatter(XPf(:,1),XPf(:,2),'red');
hold off;
