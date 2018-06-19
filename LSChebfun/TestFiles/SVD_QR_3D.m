domain = [-1 1;-1 1;-1 1];
degs = [14 14 14];
fine_degs = 2*degs;
f  = @(x,y,z) cos(3*x-sin(2*y)+exp(z));

x_f = chebpts(fine_degs(1),domain(1,:));
y_f = chebpts(fine_degs(2),domain(2,:));
z_f = chebpts(fine_degs(2),domain(3,:));

x_in = (1:degs(1))';
y_in = (1:degs(2))';
z_in = (1:degs(3))';

[X_f,Y_f,Z_f] = ndgrid(x_f,y_f,z_f);
XP_f = [X_f(:) Y_f(:) Z_f(:)];

ind_f = sqrt(X_f.^2+Y_f.^2+Z_f.^2)<=1;
XP_f = XP_f(ind_f,:);

[X_in,Y_in,Z_in] = ndgrid(x_in,y_in,z_in);

IND_N_2 = sqrt(X_in.^2+Y_in.^2+Z_in.^2);

ind_c = IND_N_2<=min(degs);

Mx_f = chebtech.clenshaw(chebpts(fine_degs(1)),eye(degs(1)));
My_f = chebtech.clenshaw(chebpts(fine_degs(2)),eye(degs(2)));
Mz_f = chebtech.clenshaw(chebpts(fine_degs(3)),eye(degs(3)));

M_f = kron(Mz_f,kron(My_f,Mx_f));
M_f = M_f(ind_f,:);


F = f(XP_f(:,1),XP_f(:,2),XP_f(:,3));

max_val = max(abs(F));

%New method where coeffs outside of sphere is set to zero
coeffs1 = zeros(prod(degs),1);
tic;
warning('off','all');
coeffs1(ind_c) = M_f(:,ind_c)\F;
warning('on','all');
toc;

norm(M_f*coeffs1-F,inf)

%old method
tic;
warning('off','all');
coeffs2 = M_f\F;
warning('on','all');
toc;

norm(M_f*coeffs2-F,inf)