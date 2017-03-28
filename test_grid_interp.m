clear;

domain = [-1 1];

x=chebpts(17,domain);
y=chebpts(33,domain);
z=chebpts(65,domain);

M = 6;
N=3;
chebpoints = cell(M,1);

chebmatrices = cell(M,2);

chebweights = cell(M,1);

deg_ind = [4 5 6];
degs = [17 33 65];

for i=1:M
    chebpoints{i} = chebpts(N);
    chebmatrices{i,1} = diffmat(N,1);
    chebmatrices{i,2} = diffmat(N,2);
    chebweights{i} = chebtech2.barywts(N);
    N = N+(N-1);
end

numb = 65;

%Simulate a splitting
xc = linspace(-1,1,numb)';
xc = xc(xc>1-0.75); 
yc = linspace(-1,1,numb)';
zc = linspace(-1,1,numb)';

grid_points = {xc,yc,zc};

[XC,YC,ZC] = ndgrid(xc,yc,zc);

XP = [XC(:) YC(:) ZC(:)];

[X,Y,Z] = ndgrid(x,y,z);

F = X.^2+Y.*X+Z.^3;

[Y2,Z2,X2] = ndgrid(y,z,x);

F2 = X2.^2+Y2.*X2+Z2.^3;

Fdx = 2*X+Y;
Fdy = X;
Fdz = 3*Z.^2;

diff_dim = 2;
H = (shiftdim(F,diff_dim-1));
perm_degs = size(H);
H = reshape(H,[perm_degs(1) prod(perm_degs(2:end))]);
H = 2*(domain(2)-domain(1))^-1*chebmatrices{deg_ind(diff_dim),1}*H;
H = reshape(H,perm_degs);
Gd = shiftdim(H,3-(diff_dim-1));


F2dx = 2*X2+Y2;
F2dy = X2;
F2dz = 3*Z2.^2;

diff_dim = 1;
H = (shiftdim(F2,diff_dim-1));
perm_degs = size(H);
H = reshape(H,[perm_degs(1) prod(perm_degs(2:end))]);
H = 2*(domain(2)-domain(1))^-1*chebmatrices{deg_ind(2),1}*H;
H = reshape(H,perm_degs);
Gd2 = shiftdim(H,3-(diff_dim-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tic;
% G = F;
% 
% h = @(x) 2/(domain(2)-domain(1))*x-(domain(2)+domain(1))/(domain(2)-domain(1));
% 
% for k=1:ndims(XC)
%     G = bary(h(grid_points{k}),G,chebpoints{deg_ind(k)},chebweights{deg_ind(k)});
%     G = permute(G,[2 3 1]);
% end
% toc
% 
% FC = XC.^2+YC.*XC+ZC.^3;
% 
% max(abs(FC(:)-G(:)))
% 
% FUNS = zeros(length(XP),1);
% 
% tic;
% for i=1:size(XP,1)
%     G = F;
%     for k=1:size(XP,2)
%         G = bary(XP(i,k),G,chebpoints{deg_ind(k)},chebweights{deg_ind(k)});
%     end
%     FUNS(i) = G;
% end
% toc
% 
% max(abs(FUNS(:)-FC(:)))




