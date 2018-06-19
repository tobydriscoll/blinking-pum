clear all;
Nx = 33;
Ny = 33;

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

Mx = cos( (0:Nx-1).*acos(XPf(:,1))); 
My = cos( (0:Nx-1).*acos(XPf(:,2))); 

%Construct Interp matrix for boundary points
tic;
M = zeros(length(XPf),length(XP));
for i=1:length(XPf)
M(i,:) = kron(My(i,:),Mx(i,:));
end
toc

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

Vc = reshape(V,Nx,Ny);

Vdiffc = [computeDerCoeffs(Vc);zeros(1,length(y))];

norm(A*Vdiffc(:)-bdx,inf)

Mx = cos( (0:Nx-1).*acos(x) ); 
My = cos( (0:Ny-1).*acos(y) ); 

M = kron(My,Mx);



V1 = reshape(M*V,Nx,Ny);



V2 = reshape(M*Vdiffc(:),Nx,Ny);

surf(X,Y,V1);
hold on;
scatter(XPf(:,1),XPf(:,2),'red');
hold off;

figure();

surf(X,Y,V2);
hold on;
scatter(XPf(:,1),XPf(:,2),'red');
hold off;

function cout = computeDerCoeffs(c)
%COMPUTEDERCOEFFS   Recurrence relation for coefficients of derivative.
%   C is the matrix of Chebyshev coefficients of a (possibly array-valued)
%   CHEBTECH object.  COUT is the matrix of coefficients for a CHEBTECH object
%   whose columns are the derivatives of those of the original.

    [n, m] = size(c);
    cout = zeros(n-1, m);                        % Initialize vector {c_r}
    w = repmat(2*(1:n-1)', 1, m);
    v = w.*c(2:end,:);                           % Temporal vector
    cout(n-1:-2:1,:) = cumsum(v(n-1:-2:1,:), 1); % Compute c_{n-2}, c_{n-4}, ...
    cout(n-2:-2:1,:) = cumsum(v(n-2:-2:1,:), 1); % Compute c_{n-3}, c_{n-5}, ...
    cout(1,:) = .5*cout(1,:);                    % Adjust the value for c_0
end