% This method sets up the linear operators to be used in the patches for
% the theta method.
%
%   input:
%    Tree: Tree that has the current time step stored
%       L: Operator of PDE u_t = L u
%       B: cell array of boundary condition operator in order NORTH SOUTH EAST WEST
%   force: forcing term
%  border: rhs of boundary condution
%  output:
%     rhs: returns RHS for the theta method.
function [rhs] = setLinOps(PUApprox,OPER,B,force,border)

rhs = zeros(length(PUApprox),1);

step = zeros(length(PUApprox.leafArray),1);

for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end

for k=1:length(PUApprox.leafArray)
    points = PUApprox.leafArray{k}.points;
    sol = force(points(:,1),points(:,2));
    dim = PUApprox.leafArray{k}.degs;
    
    %Determine the outer boundary points for each side (north south east
    %west) and the inner boundary.
    [out_border_s,~,in_border,~] = FindBoundaryIndex2DSides(dim,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox);
    
    
    Dx = kron(eye(dim(2)),diffmat(dim(1),1,PUApprox.leafArray{k}.domain(1,:)));
    Dy = kron(diffmat(dim(2),1,PUApprox.leafArray{k}.domain(2,:)),eye(dim(1)));
    
    Dxx = kron(eye(dim(2)),diffmat(dim(1),2,PUApprox.leafArray{k}.domain(1,:)));
    Dyy = kron(diffmat(dim(2),2,PUApprox.leafArray{k}.domain(2,:)),eye(dim(1)));
    
    E = eye(prod(dim));

    %Determine operator
    OP = OPER(E,points(:,1),points(:,2),Dx,Dy,Dxx,Dyy);
    
    %Determine outer boundary operators and rhs for each of the sides.
    for i=1:4
        if any(out_border_s{i}) && ~isempty(B{i})
            OP(out_border_s{i},:) = ...
                B{i}(E(out_border_s{i},:),diag(points(out_border_s{i},1)),diag(points(out_border_s{i},2)),Dx(out_border_s{i},:),Dy(out_border_s{i},:),Dxx(out_border_s{i},:),Dyy(out_border_s{i},:));
            sol(out_border_s{i}) = border{i}(points(out_border_s{i},1),points(out_border_s{i},2));
        end
    end
    
    %Enforce interface conditions
    sol(in_border) = 0;
    OP(in_border,:) = E(in_border,:);
    
    
    rhs(step(k)+(1:length(PUApprox.leafArray{k}))) = sol;

    
    PUApprox.leafArray{k}.linOp = OP;
    
    if length(PUApprox.leafArray)>1
        
        [L,U,p] = lu(OP,'vector');
        PUApprox.leafArray{k}.L = L;
        PUApprox.leafArray{k}.U = U;
        PUApprox.leafArray{k}.p = p;
        
    end
    
    if ~PUApprox.ChebRoot.is_leaf
                PUApprox.leafArray{k}.Binterp = PUApprox.ChebRoot.interpSparseMatrixZone(points(in_border,:));
    end
end

end

