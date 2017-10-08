overlap = 0.1;

domain = [-1 1;-1 1];

deg_in = [5 5];
split_flag = [1 1];
tol = 1e-6;
maxit = 1200;
cheb_length  = 33;
dim = [33 33];

Tree = ChebPatch(domain,domain,domain,deg_in,split_flag,tol);

Tree = Tree.split(1);
Tree.split();

Tree.split();
%Tree.split();

%Tree.split();

%Here split will split every leaf of the tree.
% Tree = Tree.split(0.1/4,1);
% Tree.split(0.1/4);
% 
% Tree.split(0.1/2);
% Tree.split(0.1/2);
% 
% Tree.split(0.1);

LEAVES = Tree.collectLeaves({});

force = @(x) ones(length(x),1);

boundry = @(x) zeros(length(x),1);

for k=1:length(LEAVES)
    
    dim = LEAVES{k}.degs;
    
    [out_border,in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
    
    Dxx = kron(eye(dim(2)),diffmat(dim(1),2,LEAVES{k}.domain(1,:)));
    Dyy = kron(diffmat(dim(2),2,LEAVES{k}.domain(2,:)),eye(dim(1)));
    
    E = eye(prod(dim));
    
    Lap = Dxx+Dyy;
    
    Lap(in_border,:) = E(in_border,:);
    Lap(out_border,:) = E(out_border,:);
    
    LEAVES{k}.linOp = Lap;
end

u_0 = zeros(length(Tree),1);

x = linspace(-1,1,100)';
[X,Y] = ndgrid(x);

x1 = chebpts(33,LEAVES{1}.domain(1,:));
y1 = chebpts(33,LEAVES{1}.domain(2,:));
[X1,Y1] = ndgrid(x1,y1);

[out_border1,in_border1] = FindBoundaryIndex2D(dim,LEAVES{1}.domain(),domain);

while true
    u_next = [];
    Tree.sample(u_0);
    
    F = Tree.evalfGrid({x x},1,0);
    
    surf(X,Y,F);
    
    BVALS = Tree.evalf([X1(in_border1) Y1(in_border1)],1,0);
    
    norm(F(in_border) - LEAVES{1}.values(in_border))
    
    
    %subplot(121)
    %pcolor(X1,Y1,F);
    %Tree.plotdomain;
    %axis([LEAVES{1}.domain(1,:) LEAVES{1}.domain(2,:)])
    
    %subplot(122)
    %pcolor(X1,Y1,LEAVES{1}.values)
    
    pause(0.01);
    
    for k=1:length(LEAVES)
        [out_border,in_border] = FindBoundaryIndex2D(dim,LEAVES{k}.domain(),domain);
        %F = Tree.evalfGrid(LEAVES{k}.leafGrids,1,0);
        rhs_k = zeros(length(LEAVES{k}),1);
        X_k = LEAVES{k}.points;
        rhs_k(~out_border & ~in_border) = force(X_k(~out_border & ~in_border,:));
        %rhs_k(out_border) = boundry(X_k(out_border,:));
        %rhs_k(in_border) = F(in_border);
        rhs_k(in_border)  = Tree.evalf(X_k(in_border,:),1,0);
        %rhs_k(in_border)  = Tree.evalfZone(X_k(in_border,:));
        u_next = [u_next;LEAVES{k}.linOp\rhs_k];
    end
    norm(u_0-u_next)
    
    if norm(u_0-u_next)<1e-6
        break
    end
    
    u_0 = u_next;
end