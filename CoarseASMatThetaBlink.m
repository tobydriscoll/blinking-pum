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
function [ Mat ] = CoarseASMatThetaBlink(PUApprox,B,theta,dt,t,lambda,dlambda_dt,alpha,gamma)

PUApprox.Coarsen();

xc2xs = @(xc) gamma*xc./(alpha^2-xc.^2);
yc2ys = @(t,yc) (yc+1).*(1+lambda(t))/2-1;  % map from [-1,1] to [-1,rho]

% Transformation from strip (xs,ys) to eye (xe,ye)
s2e = @(xs,ys) deal(real(tanh((xs+1i*ys)/2)),imag(tanh((xs+1i*ys)/2)));
lapfactor = @(x,y) 1./(cos(y) + cosh(x)).^2; 

c2e = @(t,xc,yc) s2e( xc2xs(xc), yc2ys(t,yc) );  % composite c->e

step = 0;

ii = [];
jj = [];
zz = [];

%y_e = @(x,y) sin(y)./(cos(y)+cosh(x));
%x_e = @(x,y) sinh(x)./(cos(y)+cosh(x));

lapfactor = @(x,y) 1./(cos(y) + cosh(x)).^2;

for k=1:length(PUApprox.leafArray)
    
    cdegs = PUApprox.leafArray{k}.cdegs;

    %Determine boundary and interface indicies
    [out_border_c,~,in_border] = FindBoundaryIndex2DSides(cdegs,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox); 
    
    pointsl = PUApprox.leafArray{k}.points();
    
    g = PUApprox.leafArray{k}.leafGrids();
    
    xc = g{1}; yc = g{2};
    
    
    index_n = (1:length(PUApprox.leafArray{k}))';
    index_n = index_n(in_border);
    
    % add them to the matrix
    % Note that ibb is the indicies of the interface border points.
    % This is why we use index_n(iib) here.
    [iib,jjb,zzb] = PUApprox.ChebRoot.interpMatrixZone_vecs(pointsl(in_border,:));
    
    ii = [ii;index_n(iib)+step];
    jj = [jj;jjb];
    zz = [zz;-zzb];    

    Dxc = diffmat(cdegs(1),1,PUApprox.leafArray{k}.domain(1,:));
    Dyc = diffmat(cdegs(2),1,PUApprox.leafArray{k}.domain(2,:));
    
    dyc_dys = @(t) 2/(1+lambda(t));
    Dys = @(t) Dyc*dyc_dys(t);
    dyc_dt = @(t) -dlambda_dt(t)*(yc+1)/(lambda(t)+1);
    
    dxs_dxc = gamma*(alpha^2+xc.^2)./(alpha.^2-xc.^2).^2;
    Dxs = diag(1./dxs_dxc)*Dxc;
    
    [xe,ye] = c2e(xc,yc,t);
    
    [Xe,Ye] = ndgrid(xe,ye);
    
        S2E = lapfactor(Xe,Ye).^(-1);
    
    Dxxe = kron(eye(cdegs(2)),S2E.*Dxs^2);
    Dyye = kron(S2E.*Dys(t+dt)^2,eye(cdegs(1)));

    OP = Dxxe+Dyye-kron(diag(dyc_dt(t+dt))*Dyc,eye(cdegs(1)));
    
    OP = eye(prod(cdegs)) - dt*theta*OP;
    
    E = eye(prod(cdegs));
    
    %Incorporate boundary conditions
    for i=1:4
        if any(out_border_c{i}) && ~isempty(B{i})
        OP(out_border_c{i},:) = ...
            B{i}(E(out_border_c{i},:),pointsl(out_border_c{i},1),pointsl(out_border_c{i},2),0,0,0,0,t+dt);
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

