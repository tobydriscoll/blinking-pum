clear all;
Nx = 35;
Ny = 35;

Nxf = 2*Nx;
Nyf = 2*Ny;

Dx = diffmat(Nx);
B_n = Nx*Nx;

tol = 1e-14;

x = chebpts(Nx);
y = chebpts(Ny);

[X,Y] = ndgrid(x,y);

x_f = chebpts(Nxf)';
y_f = chebpts(Nyf)';

[Xf,Yf] = ndgrid(x_f,y_f);

XP = [X(:) Y(:)];
XPf = [Xf(:) Yf(:)];

grid_sq_ind = (XPf(:,1)-1).^2+(XPf(:,2)-1).^2<=4;

%fine grid in domain
XPf = XPf(grid_sq_ind,:);

theta = linspace(pi,3*pi/2,B_n)';

% boundry of disk with radius 2 centered at (1,1) 
% intersected with square.
B = [2*sin(theta)+1 2*cos(theta)+1]; 
XPf = [XPf;B];

Mx = barymat(XPf(:,1),x);
My = barymat(XPf(:,2),y);


%Construct Interp matrix for boundary points
M = zeros(length(XPf),length(XP));
for i=1:length(XPf)
M(i,:) = kron(My(i,:),Mx(i,:));
end



num_in_circ = sum(grid_sq_ind);

%This is a simply matrix that
% 1) restricts the points in the grid inside
%    the domain to there function values
% 2) enforces values at the boundary
A = M;


%f = @(x,y) cos((x-1).^2+(y-1).^2);
%fdx = @(x,y) -2*(x-1).*sin((x-1).^2+(y-1).^2);

f = @(x,y) atan(x+y);
fdx = @(x,y) 1./(1+(x+y).^2);

b = f(XPf(:,1),XPf(:,2));
bdx = fdx(XPf(:,1),XPf(:,2));

tic;
V = A\b;
toc

norm(A*V-b,inf)

V1 = reshape(V,Nx,Ny);

V2 = Dx*V1;

norm(A*V2(:)-bdx,inf)

surf(X,Y,V1);
hold on;
scatter(XPf(:,1),XPf(:,2),'red');
hold off;

figure();

surf(X,Y,V2);
hold on;
scatter(XPf(:,1),XPf(:,2),'red');
hold off;