domain = [-1 1;-1 1;-1 1];
in_domain = Ball(2,[1 1 1]);
degs = [14 14 14];
mult = 2;
fine_degs = mult*degs;
c = 0.2;
bound = @(x,y,z) atan(c*(x+y+z));
f = @(x,y,z) -6*c^3*(x+y+z)./(1+c^2*(x+y+z).^2).^2;

B_n = mult*17;
b_tau = 1; %multiplier for boundary residual

% solve lap(u) = f
%           u  = bound on boundary
% 
% This is a template for how to solve BVPs in 3D with the LS square method
% 
%

x_in = (1:degs(1))';
y_in = (1:degs(2))';
z_in = (1:degs(3))';

[X_in,Y_in,Z_in] = ndgrid(x_in,y_in,z_in);

ind_c = sqrt(X_in.^2+Y_in.^2+Z_in.^2)<=max(degs);

Dx = ChebDiff(degs(1));
Dy = ChebDiff(degs(2));
Dz = ChebDiff(degs(3));

DXX = kron(eye(degs(3)),kron(eye(degs(2)),Dx^2));
DYY = kron(eye(degs(3)),kron(Dy^2,eye(degs(1))));
DZZ = kron(Dz^2,kron(eye(degs(2)),eye(degs(1))));

x = chebpts(degs(1));
y = chebpts(degs(2));
z = chebpts(degs(3));

[X,Y,Z] = ndgrid(x,y,z);

x_f = chebpts(fine_degs(1));
y_f = chebpts(fine_degs(2));
z_f = chebpts(fine_degs(3));

[Xf,Yf,Zf] = ndgrid(x_f,y_f,z_f);

XPf = [Xf(:) Yf(:) Zf(:)];

grid_sq_ind = in_domain.Interior(XPf);

%fine grid in domain
XPf = XPf(grid_sq_ind,:);

phi = linspace(pi,3*pi/2,B_n)';
theta = linspace(pi-pi/2,3*pi/2-pi/2,B_n)';

[THETA,PHI] = ndgrid(theta,phi);

% boundry of disk with radius 2 centered at (1,1) 
% intersected with square.
B = [2*sin(THETA(:)).*cos(PHI(:)) 2*sin(THETA(:)).*sin(PHI(:)) 2*cos(THETA(:))]+1; 
%XPf = [XPf;B];

Mx = chebtech.clenshaw(chebpts(degs(1)*2),eye(degs(1)));
My = chebtech.clenshaw(chebpts(degs(2)*2),eye(degs(2)));
Mz = chebtech.clenshaw(chebpts(degs(3)*2),eye(degs(3)));

M = kron(Mz,kron(My,Mx));

M = M(grid_sq_ind,:);

%Construct Interp matrix for boundary points
MB = zeros(length(B),prod(degs));

MBx = clenshaw(B(:,1),eye(degs(1)));
MBy = clenshaw(B(:,2),eye(degs(2)));
MBz = clenshaw(B(:,3),eye(degs(3)));

for i=1:length(B)
MB(i,:) = kron(MBz(i,:),kron(MBy(i,:),MBx(i,:)));
end

grid_sq_ind_b =  XPf(:,1)==1 | XPf(:,2)==1 | XPf(:,3)==1 | ((XPf(:,1)-1).^2+(XPf(:,2)-1).^2+(XPf(:,3)-1).^2==4);

grid_sq_ind_in = ~grid_sq_ind_b;
num_in_circ = sum(grid_sq_ind);

%This is a simply matrix that
% 1) restricts the points in the grid inside
%    the domain to there function values
% 2) enforces values at the boundary
A = zeros(length(XPf),prod(degs));

A(grid_sq_ind_in,:) = M(grid_sq_ind_in,:)*(DXX+DYY+DZZ);
A(grid_sq_ind_b,:) = M(grid_sq_ind_b,:);

A = [A;b_tau*MB];

b = zeros(length(XPf),1);
b(grid_sq_ind_in,:) = f(XPf(grid_sq_ind_in,1),XPf(grid_sq_ind_in,2),XPf(grid_sq_ind_in,3));
b(grid_sq_ind_b,:) = bound(XPf(grid_sq_ind_b,1),XPf(grid_sq_ind_b,2),XPf(grid_sq_ind_b,3));

b = [b;b_tau*bound(B(:,1),B(:,2),B(:,3))];


bf = bound(XPf(:,1),XPf(:,2),XPf(:,3));

tic;
V1 = A\b;
toc

V2 = zeros(prod(degs),1);
%basic solution
tic;
V2(ind_c) = A(:,ind_c)\b;
toc

norm(A*V1-b,inf)
norm(M*V1-bf,inf)

norm(A(:,ind_c)*V2(ind_c)-b,inf)
norm(M*V2-bf,inf)
