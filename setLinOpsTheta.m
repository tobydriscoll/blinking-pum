% This method sets up the linear operators to be used in the patches for
% the theta method.
%
%   input:
%    Tree: Tree that has the current time step stored
%       L: Operator of PDE u_t = L u
%       B: cell array of boundary condition operator in order NORTH SOUTH EAST WEST
%  border: rhs of boundary condution
%   theta: parameter used in theta method
%    dt,t: time step and current time
%  output:
%     rhs: returns RHS for the theta method.
function [rhs] = setLinOpsTheta(Tree,L,B,border,theta,dt,t)

if Tree.is_leaf
    LEAVES = {Tree};
else
    LEAVES = Tree.collectLeaves();    
end

rhs = zeros(length(Tree),1);

step = zeros(length(LEAVES),1);

for k=2:length(LEAVES)
    step(k) = step(k-1) + length(LEAVES{k-1});
end

for k=1:length(LEAVES)
    points = LEAVES{k}.points;

    dim = LEAVES{k}.degs;
    
    [out_border_c,~,in_border,~] = FindBoundaryIndex2DSides(dim,LEAVES{k}.domain,LEAVES{k}.outerbox);
    
    Dx = kron(eye(dim(2)),diffmat(dim(1),1,LEAVES{k}.domain(1,:)));
    Dy = kron(diffmat(dim(2),1,LEAVES{k}.domain(2,:)),eye(dim(1)));
    
    Dxx = kron(eye(dim(2)),diffmat(dim(1),2,LEAVES{k}.domain(1,:)));
    Dyy = kron(diffmat(dim(2),2,LEAVES{k}.domain(2,:)),eye(dim(1)));
    
    E = eye(prod(dim));

    OP = L(E,points(:,1),points(:,2),Dx,Dy,Dxx,Dyy,t);
    
    sol = dt*(theta-1)*OP*LEAVES{k}.values(:);
    
    OP = E - dt*theta*OP;
    
    for i=1:4
        if any(out_border_c{i}) && ~isempty(B{i})
            OP(out_border_c{i},:) = ...
                B{i}(E(out_border_c{i},:),diag(points(out_border_c{i},1)),diag(points(out_border_c{i},2)),Dx(out_border_c{i},:),Dy(out_border_c{i},:),Dxx(out_border_c{i},:),Dyy(out_border_c{i},:));
            sol(out_border_c{i}) = border{i}(points(out_border_c{i},1),points(out_border_c{i},2));
        end
    end
    
    
    sol(in_border) = 0;
    
    rhs(step(k)+(1:length(LEAVES{k}))) = sol;
    
    OP(in_border,:) = E(in_border,:);
    
    LEAVES{k}.linOp = OP;
    
    if ~Tree.is_leaf
                LEAVES{k}.Binterp = Tree.interpSparseMatrixZone(points(in_border,:));
    end
end

end

