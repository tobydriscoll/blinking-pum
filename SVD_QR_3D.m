domain = [-1 1;-1 1;-1 1];
in_domain = Ball(1,[0 0 0]);
degs = [17 17 17];
fine_degs = 2*degs;
tol = 1e-6;
dim  = 2;
f  = @(x,y,z) cos(3*x-sin(2*y)+exp(z));
tau1  = 1e-3;
%f = @(x,y) x;

%LSPATCH = LSPatch2D('InnerDomain',in_domain,'domain',domain,'degreeIndex',[5 5]);

x_f = chebpts(fine_degs(1),domain(1,:));
y_f = chebpts(fine_degs(2),domain(2,:));
z_f = chebpts(fine_degs(2),domain(3,:));

x = chebpts(degs(1),domain(1,:));
y = chebpts(degs(2),domain(2,:));
z = chebpts(degs(2),domain(3,:));

[X_f,Y_f,Z_f] = ndgrid(x_f,y_f,z_f);
[X,Y,Z] = ndgrid(x,y,z);

XP_f = [X_f(:) Y_f(:) Z_f(:)];
XP = [X(:) Y(:) Z(:)];

ind_f = in_domain.Interior(XP_f);
ind = in_domain.Interior(XP);
XP_f = XP_f(ind_f,:);
XP = XP(ind,:);

Mx_f = clenshaw(chebpts(fine_degs(1)),eye(degs(1)));
My_f = clenshaw(chebpts(fine_degs(2)),eye(degs(2)));
Mz_f = clenshaw(chebpts(fine_degs(3)),eye(degs(3)));

Mx = clenshaw(chebpts(degs(1)),eye(degs(1)));
My = clenshaw(chebpts(degs(2)),eye(degs(2)));
Mz = clenshaw(chebpts(degs(3)),eye(degs(3)));

M_f = kron(Mz_f,kron(My_f,Mx_f));
M_f = M_f(ind_f,:);

M = kron(Mz,kron(Mx,My));
M = M(ind,:);

F = f(XP_f(:,1),XP_f(:,2),XP_f(:,3));
F_l1 = f(XP(:,1),XP(:,2),XP(:,3));

max_val = max(abs(F));

tic;
warning('off','all');
coeffs1 = M_f\F;
warning('on','all');
toc;
norm(M_f*coeffs1-F,inf)


% 
% x_ind = 0:(degs(1)-1);
% y_ind = 0:(degs(2)-1);
% [X_in,Y_in] = ndgrid(x_ind,y_ind);

% subplot(2,2,1);
% surf(X_in,Y_in,reshape(abs(coeffs1),degs));
% title('basic solution');
% 
% subplot(2,2,2);
% surf(X_in,Y_in,reshape(abs(coeffs2),degs));
% title('least squares solution');
% 
% subplot(2,2,3);
% surf(X_in,Y_in,reshape(abs(coeffs3),degs));
% title('magic l_1 solution solution');