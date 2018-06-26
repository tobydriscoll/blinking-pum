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
function [rhs] = setLinOpsThetaBlink(PUfun,B,border,theta,dt,t,lambda,dlambda_dt,alpha,gamma)


xc2xs = @(xc) gamma*xc./(alpha^2-xc.^2);
yc2ys = @(t,yc) (yc+1).*(1+lambda(t))/2-1;  % map from [-1,1] to [-1,rho]

% Transformation from strip (xs,ys) to eye (xe,ye)
s2e = @(xs,ys) deal(real(tanh((xs+1i*ys)/2)),imag(tanh((xs+1i*ys)/2)));
lapfactor = @(x,y) 1./(cos(y) + cosh(x)).^2; 

c2e = @(t,xc,yc) s2e( xc2xs(xc), yc2ys(t,yc) );  % composite c->e

rhs = zeros(length(PUfun),1);

step = zeros(length(PUfun.leafArray),1);

for k=2:length(PUfun.leafArray)
    step(k) = step(k-1) + length(PUfun.leafArray{k-1});
end

for k=1:length(PUfun.leafArray)
    points = PUfun.leafArray{k}.points;

    degs = PUfun.leafArray{k}.degs;
    
    [out_border_c,~,in_border,~] = FindBoundaryIndex2DSides(degs,PUfun.leafArray{k}.domain,PUfun.leafArray{k}.outerbox);
    
    g = PUfun.leafArray{k}.leafGrids();
    
    xc = g{1}; yc = g{2};
    
    Dxc = diffmat(degs(1),1,PUfun.leafArray{k}.domain(1,:));
    Dyc = diffmat(degs(2),1,PUfun.leafArray{k}.domain(2,:));
    
    dyc_dys = @(t) 2/(1+lambda(t));
    Dys = @(t) Dyc*dyc_dys(t);
    dyc_dt = @(t) -dlambda_dt(t)*(yc+1)/(lambda(t)+1);

    dxs_dxc = gamma*(alpha^2+xc.^2)./(alpha.^2-xc.^2).^2;
    Dxs = diag(1./dxs_dxc)*Dxc;
    
    [xe,ye] = c2e(xc,yc,t);
    
    [Xe,Ye] = ndgrid(xe,ye);
    
    S2E = lapfactor(Xe,Ye).^(-1);
    
    
    Dxxe = kron(eye(degs(2)),S2E.*Dxs^2);
    Dyye = kron(S2E.*Dys(t)^2,eye(degs(1)));
    
    OP = Dxxe+Dyye+kron(diag(dyc_dt(t))*Dyc,eye(degs(1)));
    
    PUfun.leafArray{k}.linOp_f = OP;
    
    sol = PUfun.leafArray{k}.values(:)+dt*(1-theta)*OP*PUfun.leafArray{k}.values(:);
    
    [xe,ye] = c2e(xc,yc,t+dt);
    
    [Xe,Ye] = ndgrid(xe,ye);
    
    S2E = lapfactor(Xe,Ye).^(-1);
    
    Dyye = kron(S2E.*Dys(t+dt)^2,eye(degs(1)));
    
    OP = Dxxe+Dyye+kron(diag(dyc_dt(t+dt))*Dyc,eye(degs(1)));
    
    OP = eye(prod(degs)) - dt*theta*OP;
    
    E = eye(prod(degs));
    
    for i=1:4
        if any(out_border_c{i}) && ~isempty(B{i})
            OP(out_border_c{i},:) = ...
                B{i}(E(out_border_c{i},:),points(out_border_c{i},1),points(out_border_c{i},2),0,0,0,0,t+dt);
            sol(out_border_c{i}) = border{i}(points(out_border_c{i},1),points(out_border_c{i},2),t+dt);
        end
    end
    
    
    sol(in_border) = 0;
    
    rhs(step(k)+(1:length(PUfun.leafArray{k}))) = sol;
    
    OP(in_border,:) = E(in_border,:);
    
    PUfun.leafArray{k}.linOp = OP;
    
    if ~PUfun.ChebRoot.is_leaf
                PUfun.leafArray{k}.Binterp = PUfun.ChebRoot.interpSparseMatrixZone(points(in_border,:));
    end
end

end

