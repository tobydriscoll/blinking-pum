clear;

domain = [-1 1];

standard_degs = [3 5 9 17 33 65];

deg_ind = [5 5 5 5];
degs = standard_degs(deg_ind);

x=chebpts(degs(1),domain);
y=chebpts(degs(2),domain);
z=chebpts(degs(3),domain);
w=chebpts(degs(4),domain);


M = 6;
N=3;
chebpoints = cell(M,1);

chebmatrices = cell(M,2);

chebweights = cell(M,1);



for i=1:M
    chebpoints{i} = chebpts(N);
    chebmatrices{i,1} = diffmat(N,1);
    chebmatrices{i,2} = diffmat(N,2);
    chebweights{i} = chebtech2.barywts(N);
    N = N+(N-1);
end

numb = 65;

%Simulate a splitting
xc = linspace(-1,1,65)';
yc = linspace(-1,1,65)';
yc = yc(yc>1-0.75); 
zc = linspace(-1,1,65)';
wc = linspace(-1,1,65)';

grid_points = {xc,yc,zc,wc};

[X2C,Y2C] = ndgrid(xc,yc);
[X3C,Y3C,Z3C] = ndgrid(xc,yc,zc);
[X4C,Y4C,Z4C,W4C] = ndgrid(xc,yc,zc,wc);

XP2 = [X2C(:) Y2C(:)];
XP3 = [X3C(:) Y3C(:) Z3C(:)];
XP4 = [X4C(:) Y4C(:) Z4C(:) W4C(:)];


[X2,Y2] = ndgrid(x,y);
[X3,Y3,Z3] = ndgrid(x,y,z);
[X4,Y4,Z4,W4] = ndgrid(x,y,z,w);

F2 = X2.^2+Y2.*X2;
F3 = X3.^2+Y3.*X3+Z3.^3;
F4 = X4.^2+Y4.*X4+Z4.^3+W4.^4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %2D times
% tic;
% G = F2;
% 
% h = @(x) 2/(domain(2)-domain(1))*x-(domain(2)+domain(1))/(domain(2)-domain(1));
% 
% for k=1:ndims(X2C)
%     G = bary(h(grid_points{k}),G,chebpoints{deg_ind(k)},chebweights{deg_ind(k)});
%     G = shiftdim(G,1);
% end
% toc
% 
% F2C = X2C.^2+Y2C.*X2C;
% 
% max(abs(F2C(:)-G(:)))
% 
% FUNS = zeros(length(XP2),1);
% 
% tic;
% for i=1:size(XP2,1)
%     G = F2;
%     for k=1:size(XP2,2)
%         G = bary(XP2(i,k),G,chebpoints{deg_ind(k)},chebweights{deg_ind(k)});
%     end
%     FUNS(i) = G;
% end
% toc
% 
% max(abs(FUNS(:)-F2C(:)))
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %3D times
tic;
G = F3;

h = @(x) 2/(domain(2)-domain(1))*x-(domain(2)+domain(1))/(domain(2)-domain(1));

for k=1:ndims(X3C)
    G = bary(h(grid_points{k}),G,chebpoints{deg_ind(k)},chebweights{deg_ind(k)});
    G = shiftdim(G,1);
end
toc

F3C = X3C.^2+Y3C.*X3C+Z3C.^3;

max(abs(F3C(:)-G(:)))

FUNS = zeros(length(XP3),1);

tic;
for i=1:size(XP3,1)
    G = F3;
    for k=1:size(XP3,2)
        G = bary(XP3(i,k),G,chebpoints{deg_ind(k)},chebweights{deg_ind(k)});
    end
    FUNS(i) = G;
end
toc

max(abs(FUNS(:)-F3C(:)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4D times

% tic;
% G = F4;
% 
% h = @(x) 2/(domain(2)-domain(1))*x-(domain(2)+domain(1))/(domain(2)-domain(1));
% 
% for k=1:ndims(X4C)
%     G = bary(h(grid_points{k}),G,chebpoints{deg_ind(k)},chebweights{deg_ind(k)});
%     G = shiftdim(G,1);
% end
% toc
% 
% F4C = X4C.^2+Y4C.*X4C+Z4C.^3+W4C.^4;
% 
% max(abs(F4C(:)-G(:)))
% 
% FUNS = zeros(length(XP3),1);
% 
% tic;
% for i=1:size(XP4,1)
%     G = F4;
%     for k=1:size(XP4,2)
%         G = bary(XP4(i,k),G,chebpoints{deg_ind(k)},chebweights{deg_ind(k)});
%     end
%     FUNS(i) = G;
% end
% toc
% 
% max(abs(FUNS(:)-F3C(:)))
