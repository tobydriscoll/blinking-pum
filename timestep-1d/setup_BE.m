function [data,DAEfun,jacfun] = setup_BE(ODE,BC,dt,domain)

Nd = length(domain);
xlim = cat(1,domain.xlim);
n = cat(1,domain.n);
idx = false(sum(n),1);

data = struct;

for d = 1:Nd
    data(d).x = chebpts(domain(d).n,domain(d).xlim);
    data(d).D = diffmat(domain(d).n,domain(d).xlim);
    D1 = data(d).D;  D2 = D1^2;
    
    res = @(u) ODE.fun(u,D1,D2);   data(d).resfun = res;
    jac = @(u) ODE.jac(u,D1,D2);   data(d).jacfun = jac;
    
    M = eye(domain(d).n); M([1 end],:) = 0;
    data(d).M = M;
    
    data(d).idx = idx;
    data(d).idx( sum(n(1:d-1))+(1:n(d)) ) = true;
    
end

% cross-interpolation matrices
for d = 1:Nd
    
    Bl = cell(1,Nd);  Br = cell(1,Nd);
    for j = 1:Nd
        if xlim(d,1) > xlim(j,1) && xlim(d,1) < xlim(j,2)
            Bl{j} = barymat(xlim(d,1),data(j).x);
        else
            Bl{j} = sparse(1,n(j));
        end
        
        if xlim(d,2) > xlim(j,1) && xlim(d,2) < xlim(j,2)
            Br{j} = barymat(xlim(d,2),data(j).x);
        else
            Br{j} = sparse(1,n(j));
        end
    end
    data(d).Bl = cell2mat(Bl);  data(d).Br = cell2mat(Br);
    
end

% subdomain functions
for d = 1:Nd
    lbc = [];  if d==1, lbc = BC(1); end
    rbc = [];  if d==Nd, rbc = BC(2); end
    Bl = data(d).Bl;
    Br = data(d).Br;
end

    function [res,jac] = DAE(u,un)
        res = cell(1,Nd);  jac = cell(1,Nd);
        for i = 1:Nd
            [res{i},jac{i}] = DAEfuns(u,un,data(i).resfun,data(i).jacfun,dt,lbc,rbc,...
                data(i).idx,data(i).M,Bl,Br);
        end
    end

DAEfun = @DAE;
jacfun = @(z,jac) jacapply(z,domain,data,jac);

end

function [resfun,jacfun] = DAEfuns(u,un,ODEfun,ODEjac,dt,lbc,rbc,idx,M,Bl,Br)

%maker = @funmaker;

%    function [resfun,jacfun] = funmaker(u,un)
        uloc = u(idx);
        resfun = @fun;
        jacfun = @jac;
        
        
        function r = fun(z)
            r = ODEfun(z+uloc);
            
            if ~isempty(lbc)
                r(1) = z(1) + uloc(1) - lbc;
            else
                r(1) = z(1) + uloc(1) - Bl*u;
            end
            
            if ~isempty(rbc)
                r(end) = z(end) + uloc(end) - rbc;
            else
                r(end) = z(end) + uloc(end) - Br*u;
            end
            
            r = M*(z+uloc-un(idx)) - dt*r;
            
        end
        
        function J = jac(z)
            J = ODEjac(z+uloc);
            J(1,:) = 0;  J(1,1) = 1;
            J(end,:) = 0;  J(end,end) = 1;
            J = M - dt*J;
        end
        
%    end

end

function Jfun = jacapply(z,domain,data,jac)

Nd = length(domain);
Jfun = @apply;

    function Jv = apply(v)
        Jv = cell(Nd,1);
        for d=1:Nd
            idx = data(d).idx;
            J = jac{d}(z(idx));
            r = zeros(domain(d).n,1);
            r(1) = data(d).Bl*v;
            r(end) = data(d).Br*v;
            Jv{d} = (J\r) - v(idx);
        end
        Jv = cell2mat(Jv);
    end
end
