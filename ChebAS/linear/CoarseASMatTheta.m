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
function [ Mat ] = CoarseASMatTheta(PUApprox,L,B,theta,dt,t)

PUApprox.Coarsen();

step = 0;

ii = [];
jj = [];
zz = [];

for k=1:length(PUApprox.leafArray)
    
    cdegs = PUApprox.leafArray{k}.cdegs;

    %Determine boundary and interface indicies
    [out_border_c,~,in_border] = FindBoundaryIndex2DSides(cdegs,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox); 
    
    pointsl = PUApprox.leafArray{k}.points();
    
    index_n = (1:length(PUApprox.leafArray{k}))';
    index_n = index_n(in_border);
    
    % add them to the matrix
    % Note that ibb is the indicies of the interface border points.
    % This is why we use index_n(iib) here.
    [iib,jjb,zzb] = PUApprox.ChebRoot.interpMatrixZone_vecs(pointsl(in_border,:));
    
    ii = [ii;index_n(iib)+step];
    jj = [jj;jjb];
    zz = [zz;-zzb];    

    Dx = kron(eye(cdegs(2)),diffmat(cdegs(1),1,PUApprox.leafArray{k}.domain(1,:)));
    Dy = kron(diffmat(cdegs(2),1,PUApprox.leafArray{k}.domain(2,:)),eye(cdegs(1)));
    
    Dxx = kron(eye(cdegs(2)),diffmat(cdegs(1),2,PUApprox.leafArray{k}.domain(1,:)));
    Dyy = kron(diffmat(cdegs(2),2,PUApprox.leafArray{k}.domain(2,:)),eye(cdegs(1)));
    
    E = eye(prod(cdegs));
    
    %Determine Operator used for the theta method
    OP = E - dt*theta*L(E,diag(pointsl(:,1)),diag(pointsl(:,2)),Dx,Dy,Dxx,Dyy,t+dt);
    
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

Mat = sparse(ii,jj,zz,length(PUApprox),length(PUApprox));

PUApprox.Refine();

end

