clear all;
Nx = 33;
Ny = 33;

Nxf = 2*Nx;
Nyf = 2*Ny;

Dx = diffmat(Nx);
Dy = diffmat(Ny);

DX = kron(eye(Ny),Dx);
DY = kron(Dy,eye(Nx));

DXX = kron(eye(Ny),Dx^2);
DYY = kron(Dy^2,eye(Nx));

B_n = Nx*Nx;

nu = 1;

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

fine_ind_in = length(XPf);

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

grid_sq_ind_b =  XPf(:,1)==1 | XPf(:,2)==1 |  ((XPf(:,1)-1).^2+(XPf(:,2)-1).^2==4);

grid_sq_ind_in = ~grid_sq_ind_b;
num_in_circ = sum(grid_sq_ind);

x0 = zeros(length(XP),1);
options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true);
fun = @(V) non_lin_pois_f(V,M,DX,DY,DXX,DYY,XPf,nu,grid_sq_ind_in,grid_sq_ind_b);
tic;
V = lsqnonlin(fun,x0,[],[],options);
toc

norm(fun(V))

Z = zeros(Nxf*Nyf,1);

Z(grid_sq_ind) = M(1:fine_ind_in,:)*V;

Z(~grid_sq_ind) = nan;

V1 = reshape(Z,Nxf,Nyf);

surf(Xf,Yf,V1);
hold on;
scatter(XPf(:,1),XPf(:,2),'red');
hold off;
