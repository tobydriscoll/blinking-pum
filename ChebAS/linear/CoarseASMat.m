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
function [ Mat_struct ] = CoarseASMat(PUApprox,Linop,B)

PUApprox.Coarsen();

step = 0;

ii = [];
jj = [];
zz = [];

for k=1:length(PUApprox.leafArray)
    cdim = PUApprox.leafArray{k}.cdegs;

    %Figure out indicies of boundary and interface points
    [out_border_s,~,in_border,~]  = FindBoundaryIndex2DSides(cdim,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox); 
    
    pointsl = PUApprox.leafArray{k}.points();
    
    index_n = (1:length(PUApprox.leafArray{k}))';
    index_n = index_n(in_border);
    
    %Figure out sparse indicies of interpolation matrix
    [iib,jjb,zzb] = PUApprox.ChebRoot.interpMatrixZone_vecs(pointsl(in_border,:));
    
    % add them to the matrix
    % Note that ibb is the indicies of the interface border points.
    % This is why we use index_n(iib) here.
    ii = [ii;index_n(iib)+step];
    jj = [jj;jjb];
    zz = [zz;-zzb];    

    Dx = kron(eye(cdim(2)),diffmat(cdim(1),1,PUApprox.leafArray{k}.domain(1,:)));
    Dy = kron(diffmat(cdim(2),1,PUApprox.leafArray{k}.domain(2,:)),eye(cdim(1)));
    
    Dxx = kron(eye(cdim(2)),diffmat(cdim(1),2,PUApprox.leafArray{k}.domain(1,:)));
    Dyy = kron(diffmat(cdim(2),2,PUApprox.leafArray{k}.domain(2,:)),eye(cdim(1)));
    
    E = eye(prod(cdim));
    
    %Determine the operator
    OP = Linop(E,pointsl(:,1),pointsl(:,2),Dx,Dy,Dxx,Dyy);
    
    %Incorporate outer boundary conditions
    for i=1:4
        if any(out_border_s{i}) && ~isempty(B{i})
        OP(out_border_s{i},:) = ...
            B{i}(E(out_border_s{i},:),pointsl(out_border_s{i},1),pointsl(out_border_s{i},2),Dx(out_border_s{i},:),Dy(out_border_s{i},:),Dxx(out_border_s{i},:),Dyy(out_border_s{i},:));
        end
    end
    
    %Incorporate interface boundary condition
    OP(in_border,:) = E(in_border,:);
    
    %Add the operator to sparse matrix
    [iid,jjd,zzd] = find(OP);
    
    ii = [ii;iid+step];
    jj = [jj;jjd+step];
    zz = [zz;zzd];
    
    step = step+prod(cdim);
end

Mat = sparse(ii,jj,zz,length(PUApprox),length(PUApprox));

[L,U,P,Q,R] = lu(Mat);

Mat_struct.L = L;
Mat_struct.U = U;
Mat_struct.P = P;
Mat_struct.Q = Q;
Mat_struct.R = R;

PUApprox.Refine();

end

