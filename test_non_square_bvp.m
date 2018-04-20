domain = [-1 1;-1 1];
in_domain = Disk(2,[1 1]);
degs = [17 17];

mult = 12;
fine_degs = mult*degs;

B_n = mult*17;

tau1 = 1e-5;
b_tau = 1;

x_in = (1:degs(1))';
y_in = (1:degs(2))';

[X_in,Y_in] = ndgrid(x_in,y_in);

ind_c = sqrt(X_in.^2+Y_in.^2)<=max(degs);

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

x_f = chebpts(fine_degs(1));
y_f = chebpts(fine_degs(2));

[Xf,Yf] = ndgrid(x_f,y_f);

XPf = [Xf(:) Yf(:)];

grid_sq_ind = (XPf(:,1)-1).^2+(XPf(:,2)-1).^2<=4;

%fine grid in domain
XPf = XPf(grid_sq_ind,:);

theta = linspace(pi,3*pi/2,B_n)';

% boundry of disk with radius 2 centered at (1,1) 
% intersected with square.
B = [2*sin(theta)+1 2*cos(theta)+1]; 
%XPf = [XPf;B];

Mx = clenshaw(x_f,eye(degs(1)));
My = clenshaw(y_f,eye(degs(2)));

M = kron(My,Mx);

M = M(grid_sq_ind,:);

%Construct Interp matrix for boundary points
MB = zeros(length(B),prod(degs));

MBx = clenshaw(B(:,1),eye(degs(1)));
MBy = clenshaw(B(:,2),eye(degs(2)));

for i=1:length(B)
MB(i,:) = kron(MBy(i,:),MBx(i,:));
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

A = [A;b_tau*MB];

f = @(x,y) zeros(size(x));
bound = @(x,y) log((x+1).^2+(y+1).^2);
fdx = @(x,y) 2*(1+x)./((1+x).^2+(1+y).^2);

b = zeros(length(XPf),1);
b(grid_sq_ind_in,:) = f(XPf(grid_sq_ind_in,1),XPf(grid_sq_ind_in,2));
b(grid_sq_ind_b,:) = bound(XPf(grid_sq_ind_b,1),XPf(grid_sq_ind_b,2));

b = [b;b_tau*bound(B(:,1),B(:,2))];


bf = bound(XPf(:,1),XPf(:,2));
bfdx = fdx(XPf(:,1),XPf(:,2));

lap = DXX+DYY+2*DXY;

tic;
V = A\b;
toc

V = zeros(prod(degs),1);
%basic solution
tic;
V(ind_c) = A(:,ind_c)\b;
toc



norm(A*V-b,inf)
norm(M*V-bf,inf)
norm(M*DX*V-bfdx,inf)

V1 = chebfun2.coeffs2vals(reshape(V,degs));

surf(X,Y,V1);
hold on;
scatter(XPf(:,1),XPf(:,2),'red');
hold on;
scatter(B(:,1),B(:,2),'blue');
hold off;
