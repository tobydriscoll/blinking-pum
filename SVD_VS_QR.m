domain = [-1 1;-1 1];
%in_domain = Disk(1,[0 0]);
in_domain = Astroid(0);
degs = [33 33];
fine_degs = 4*degs;
tol = 1e-12;
dim  = 2;
c = 1.2;
%f = @(x,y)log(c-(x.^2+y.^2));
f = @(x,y) atan(x)+y;
tau1 = 1e-7;
%f = @(x,y) x;

D_x = ChebDiff(degs(1));
D_y = ChebDiff(degs(2));
D_xx = ChebDiff(degs(1))^2;
D_yy = ChebDiff(degs(2))^2;

Lap = kron(D_yy,eye(degs(1)))+kron(eye(degs(2)),D_xx)+2*kron(D_y,D_x);

%Lap = kron(eye(degs(2)),eye(degs(1)));

%LSPATCH = LSPatch2D('InnerDomain',in_domain,'domain',domain,'degreeIndex',[4 4]);

x_f = chebpts(fine_degs(1),domain(1,:));
y_f = chebpts(fine_degs(2),domain(2,:));

x = chebpts(degs(1),domain(1,:));
y = chebpts(degs(2),domain(2,:));

x_r = domain(1,1) + diff(domain(1,:)).*rand(100,1);
y_y = domain(2,1) + diff(domain(2,:)).*rand(100,1);
%Lap = kron(diag(w_y),diag(w_x));

[X_f,Y_f] = ndgrid(x_f,y_f);
[X,Y] = ndgrid(x,y);

XP_f = [X_f(:) Y_f(:)];
XP = [X(:),Y(:)];

ind_f = in_domain.Interior(XP_f);
ind = in_domain.Interior(XP);
XP_f = XP_f(ind_f,:);
XP = XP(ind,:);

Mx_f = chebtech.clenshaw(chebpts(fine_degs(1)),eye(degs(1)));
My_f = chebtech.clenshaw(chebpts(fine_degs(2)),eye(degs(2)));

Mx = chebtech.clenshaw(chebpts(degs(1)),eye(degs(1)));
My = chebtech.clenshaw(chebpts(degs(2)),eye(degs(2)));

M_f = kron(My_f,Mx_f);
M_f = M_f(ind_f,:);

M = kron(My,Mx);
M = M(ind,:);

F = f(XP_f(:,1),XP_f(:,2));
%F_lap = lapf(XP_f(:,1),XP_f(:,2));
F_l1 = f(XP(:,1),XP(:,2));

max_val = max(abs(F));

%basic solution
tic;
coeffs1 = M_f\F;
toc;
norm(M_f*coeffs1-F,inf)/norm(F,inf)
%norm(M_f*Lap*coeffs1-F_lap,inf)/norm(F_lap,inf)

%min l_2 solution
tic;
coeffs2 = ResitrictedLS(M_f,F,eye(prod(degs)));
toc;
norm(M_f*coeffs2-F,inf)/norm(F,inf)
%norm(M_f*Lap*coeffs2-F_lap,inf)/norm(F_lap,inf)

A = partialDCT(ind_f,degs,fine_degs);
At = partialDCTtt(ind_f,degs,fine_degs);


%Tikhonov regularization
M_f2 = [tau1*Lap;M_f];
F2 = [zeros(length(coeffs1),1);F];

tic;
coeffs3 = M_f2\F2;
toc;

norm(M_f*coeffs3-F,inf)/norm(F,inf)
%norm(M_f*Lap*coeffs3-F_lap,inf)/norm(F_lap,inf)

x_ind = 0:(degs(1)-1);
y_ind = 0:(degs(2)-1);
[X_in,Y_in] = ndgrid(x_ind,y_ind);

subplot(2,2,1);
surf(X_in,Y_in,log10(reshape(abs(coeffs1),degs)));
title('basic solution');

subplot(2,2,2);
surf(X_in,Y_in,log10(reshape(abs(coeffs2),degs)));
title('min l_2 solution');

subplot(2,2,3);
surf(X_in,Y_in,log10(reshape(abs(coeffs3),degs)));
title('Tikhonov regularization');

CC = chebfun2.vals2coeffs(f(X,Y));

subplot(2,2,4);
surf(X_in,Y_in,log10(abs(CC)));
title('True Coeffs');

C1 = Chop(reshape(coeffs1,degs),tol,max_val);
C2 = Chop(reshape(coeffs2,degs),tol,max_val);
C3 = Chop(reshape(coeffs3,degs),tol,max_val);
C4 = Chop(CC,tol,max_val);

C1
C2
C3
C4

% figure();
% 
% subplot(2,2,1);
% surf(X,Y,chebfun2.coeffs2vals(reshape(coeffs1,degs)));
% title('basic solution');
% 
% subplot(2,2,2);
% surf(X,Y,chebfun2.coeffs2vals(reshape(coeffs2,degs)));
% title('min l_2 solution');
% 
% subplot(2,2,3);
% surf(X,Y,chebfun2.coeffs2vals(reshape(coeffs3,degs)));
% title('Tikhonov regularization');
% 
% CC = chebfun2.vals2coeffs(real(f(X,Y)));
% 
% subplot(2,2,4);
% surf(X,Y,f(X,Y));
% title('True SOL');

function cutoff = Chop(coeffs,tol,MaxV)

loc_tol = tol^(7/8);

for k=1:2
    
    Vals = chebfun2.coeffs2vals(coeffs);
    colChebtech = chebfun3t.unfold(Vals, k);
    fCol = chebtech2(colChebtech);
    
    
    cutoff(k) = length(simplify(fCol, tol))+1;
    
    
end

end
