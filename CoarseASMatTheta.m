% This method sets up the linear operators to be used in the patches for
% the theta method.
%
%   input:
%    Tree: Tree that has the current time step stored
%       L: Operator of PDE u_t = L u
%       B: cell array of boundary condition operator in order NORTH SOUTH EAST WEST
%   theta: parameter used in theta method
%    dt,t: time step and current time
%  output:
%     MAT: returns the coarse matrix for the theta method
function [ Mat ] = CoarseASMatTheta( Tree,L,B,theta,dt,t)

if Tree.is_leaf
    LEAVES = {Tree};
else
    LEAVES = Tree.collectLeaves();    
end

Tree.Coarsen();

step = 0;

ii = [];
jj = [];
zz = [];

for k=1:length(LEAVES)
    
    cdegs = LEAVES{k}.cdegs;

    %Determine boundary and interface indicies
    [out_border_c,~,in_border] = FindBoundaryIndex2DSides(cdegs,LEAVES{k}.domain,LEAVES{k}.outerbox); 
    
    pointsl = LEAVES{k}.points();
    
    index_n = (1:length(LEAVES{k}))';
    index_n = index_n(in_border);
    
    % add them to the matrix
    % Note that ibb is the indicies of the interface border points.
    % This is why we use index_n(iib) here.
    [iib,jjb,zzb] = Tree.interpMatrixZone_vecs(pointsl(in_border,:));
    
    ii = [ii;index_n(iib)+step];
    jj = [jj;jjb];
    zz = [zz;-zzb];    

    Dx = kron(eye(cdegs(2)),diffmat(cdegs(1),1,LEAVES{k}.domain(1,:)));
    Dy = kron(diffmat(cdegs(2),1,LEAVES{k}.domain(2,:)),eye(cdegs(1)));
    
    Dxx = kron(eye(cdegs(2)),diffmat(cdegs(1),2,LEAVES{k}.domain(1,:)));
    Dyy = kron(diffmat(cdegs(2),2,LEAVES{k}.domain(2,:)),eye(cdegs(1)));
    
    E = eye(prod(cdegs));
    
    %Determine Operator used for the theta method
    OP = E - dt*theta*L(E,pointsl(:,1),pointsl(:,2),Dx,Dy,Dxx,Dyy,t);
    
    %Incorporate boundary conditions
    for i=1:4
        if any(out_border_c{i}) && ~isempty(B{i})
        OP(out_border_c{i},:) = ...
            B{i}(E(out_border_c{i},:),pointsl(out_border_c{i},1),pointsl(out_border_c{i},2),Dx(out_border_c{i},:),Dy(out_border_c{i},:),Dxx(out_border_c{i},:),Dyy(out_border_c{i},:),t);
        end
    end
    
    %Enforce boundarty condition at interface
    OP(in_border,:) = E(in_border,:);
    
    [iid,jjd,zzd] = find(OP);
    
    ii = [ii;iid+step];
    jj = [jj;jjd+step];
    zz = [zz;zzd];
    
    step = step+prod(cdegs);
end

Mat = sparse(ii,jj,zz,length(Tree),length(Tree));

Tree.Refine();

end

