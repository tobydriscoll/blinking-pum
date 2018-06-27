function [domain,coarse,BEsys] = setup(ODE,domain)

Nd = length(domain);
xlim = cat(1,domain.xlim);
n = cat(1,domain.n);
idx = false(sum(n),1);


% This is the generic definition of the nonlinear system that defines a
% Backward Euler time step.
BEsys = @BEstep;
function [r,J] = BEstep(u,ODEfun,ODEjac,BC,un,dt)
    m = length(u);
    r = (u-un) - dt*ODEfun(u);

    r(1) = u(1) - BC(1);
    r(m) = u(m) - BC(2);
    
    J = spdiags([0;ones(m-2,1);0],0,m,m) - dt*ODEjac(u);
    J(1,:) = 0;  J(1,1) = 1;
    J(end,:) = 0;  J(end,end) = 1;
end

% points, diffmats, and subdomain problems
for d = 1:Nd
    domain(d).x = chebpts(domain(d).n,domain(d).xlim);
    domain(d).D = diffmat(domain(d).n,domain(d).xlim);
    D1 = domain(d).D;  D2 = D1^2;
    
    res = @(u) ODE.fun(u,D1,D2);   domain(d).resfun = res;
    jac = @(u) ODE.jac(u,D1,D2);   domain(d).jacfun = jac;
    
    M = speye(domain(d).n); M([1 end],:) = 0;
    domain(d).M = M;
    
    domain(d).idx = idx;
    domain(d).idx( sum(n(1:d-1))+(1:n(d)) ) = true;
    
end

% cross-interpolation matrices
for d = 1:Nd
    
    Bl = cell(1,Nd);  Br = cell(1,Nd);
    for j = 1:Nd
        if xlim(d,1) > xlim(j,1) && xlim(d,1) < xlim(j,2)
            Bl{j} = barymat(xlim(d,1),domain(j).x);
        else
            Bl{j} = sparse(1,n(j));
        end
        
        if xlim(d,2) > xlim(j,1) && xlim(d,2) < xlim(j,2)
            Br{j} = barymat(xlim(d,2),domain(j).x);
        else
            Br{j} = sparse(1,n(j));
        end
    end
    domain(d).Bl = cell2mat(Bl);  domain(d).Br = cell2mat(Br);
    
end

% partition of unity
shape = @(x) min(realmax,exp(-1./(1-x.^2))).*(abs(x)<1);

for d = 1:Nd
    if d==1
        domain(d).bump = @(x) shape((x+1)/(xlim(d,2)+1));
    elseif d==Nd
        domain(d).bump = @(x) shape((x-xlim(d,1))/(1-xlim(d,1))-1);
    else
        c = mean(xlim(d,:));
        domain(d).bump = @(x) shape((x-c)/(xlim(d,2)-c));
    end
end

function s = bumpsum(x)
    s = 0;
    for dd = 1:Nd
        s = s + domain(dd).bump(x);
    end
end
%denom = @(x) sum( cellfun(@(f) feval(f,x),{data.bump}) );

for d = 1:Nd
    domain(d).weight = @(x) domain(d).bump(x)./bumpsum(x);
end

% coarsening operations (restrict and prolong)
nc = min(400,ceil(sum(n)*0.2));
coarse.n = nc;
xc = chebpts(nc);
coarse.x = xc;
coarse.D = diffmat(nc);
coarse.resfun = @(u) ODE.fun(u,coarse.D,coarse.D^2);   
coarse.jacfun = @(u) ODE.jac(u,coarse.D,coarse.D^2);
coarse.M = spdiags([0;ones(nc-2,1);0],0,nc,nc);

for d = 1:Nd
    idx = (xc >= domain(d).xlim(1)) & (xc <= domain(d).xlim(2));
    R = barymat(xc(idx),domain(d).x);
    R = R.*domain(d).weight(domain(d).x).';
    domain(d).R = R;
    domain(d).Ridx = idx;
    
    domain(d).P = barymat(domain(d).x,xc);
end

coarse.restrict = @restrict;
    function uc = restrict(u)
        uc = zeros(coarse.n,1);
        for d = 1:Nd
            idx = domain(d).Ridx;
            uc(idx) = uc(idx) + domain(d).R*u(domain(d).idx);
        end
    end

coarse.prolong = @prolong;
    function u = prolong(uc)
        u = zeros(sum(n),1);
        for d = 1:Nd
            u(domain(d).idx) = domain(d).P*uc;
        end
    end

end