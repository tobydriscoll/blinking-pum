function [data,DAEfun,jacfun] = setup_BE(ODE,domain)

Nd = length(domain);
xlim = cat(1,domain.xlim);
n = cat(1,domain.n);
idx = false(sum(n),1);

data = struct;

% points, diffmats, and subdomain problems
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

end