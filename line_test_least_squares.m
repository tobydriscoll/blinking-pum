clear all;

Nx = 33;
Ny = 33;

Nxf = 2*Nx;
Nyf = 2*Ny;

Nxf2 = 4*Nx;
Nyf2 = 4*Ny;

Dx = diffmat(Nx);
B_n = Nx*Ny;

tol = 1e-14;

C_s = linspace(0,2);

%C_s = 0;

b = 1;
f = @(x,y) cos((x-1).^2+(y-1).^2);
fdx = @(x,y) -2*(x-1).*sin((x-1).^2+(y-1).^2);

ERR = zeros(size(C_s));
DIFF_ERR = zeros(size(C_s));

for j=1:length(C_s)
    
    A = 1; B = 1; C = C_s(j);
    LINE = @(x,y) A*x+B*y+C;
    
    LINE2 = @(x) -(C+A*x)/B;
    

    
    x = linspace(-1,1,Nx)';
    y = linspace(-1,1,Ny)';
    
    [X,Y] = ndgrid(x,y);
    
    x_f = linspace(-1,1,Nxf)';
    y_f = linspace(-1,1,Nyf)';
    
    %x_f = chebpts(Nxf)';
    %y_f = chebpts(Nyf)';
    
    [Xf,Yf] = ndgrid(x_f,y_f);
    
    x_f2 = linspace(-1,1,Nxf2)';
    y_f2 = linspace(-1,1,Nyf2)';
    
    [Xf2,Yf2] = ndgrid(x_f2,y_f2);
    
    XP = [X(:) Y(:)];
    XPf = [Xf(:) Yf(:)];
    XPf2 = [Xf2(:) Yf2(:)];
    
    grid_sq_ind = LINE(XPf(:,1),XPf(:,2))*sign(B)>=0;
    grid_sq_ind2 = LINE(XPf2(:,1),XPf2(:,2))*sign(B)>=0;
    
    %fine grid in domain
    XPf = XPf(grid_sq_ind,:);
    
   % theta = linspace(-1,(B-C)/A,2*Nx)';
   % BOUND = [theta LINE2(theta)];
   % XPf = [XPf;BOUND];
    
    XPf2 = XPf2(grid_sq_ind2,:);
    
    theta = linspace(pi,3*pi/2,B_n)';
    
    Mx = cos( (0:Nx-1).*acos(XPf(:,1)));
    My = cos( (0:Ny-1).*acos(XPf(:,2)));
    
    %Construct Interp matrix for boundary points
    M = zeros(length(XPf),length(XP));
    for i=1:length(XPf)
        M(i,:) = kron(My(i,:),Mx(i,:));
    end
    
    Mx2 = cos( (0:Nx-1).*acos(XPf2(:,1)));
    My2 = cos( (0:Ny-1).*acos(XPf2(:,2)));
    
    %Construct Interp matrix for boundary points
    M2 = zeros(length(XPf2),length(XP));
    for i=1:length(XPf2)
        M2(i,:) = kron(My2(i,:),Mx2(i,:));
    end
    
    num_in_circ = sum(grid_sq_ind);
    
    b = f(XPf(:,1),XPf(:,2));
    b2 = f(XPf2(:,1),XPf2(:,2));
    bdx = fdx(XPf2(:,1),XPf2(:,2));
    
    AM = M;
    bM = b;
    
    V = (AM)\(bM);
    
    ERR(j) = norm(M2*V-b2,inf);
    
    Vc = reshape(V,Nx,Ny);
    
    Vdiffc = [computeDerCoeffs(Vc);zeros(1,length(y))];
    
    DIFF_ERR(j) = norm(M2*Vdiffc(:)-bdx,inf);
    
end

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