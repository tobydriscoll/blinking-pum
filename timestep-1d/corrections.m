function [z,jacfun] = corrections(u,un,dt,BC,domain,data)

Nd = length(domain);
n = cat(1,domain.n);

z = cell(Nd,1);
opt = [20,-1,.5,0];

corres = cell(1,Nd);  corjac = cell(1,Nd);
for i = 1:Nd
    if i==1, lbc = BC(1); else lbc = []; end
    if i==Nd, rbc = BC(2); else rbc = []; end
    [corres{i},corjac{i}] = DAEfuns(u,un,data(i),dt,lbc,rbc);
end

for i = 1:Nd
    f = nsoldbridge(corres{i},corjac{i});
    [z{i},ithist] = nsold(zeros(n(i),1),f,[1e-8 1e-8],opt);
end
z = cell2mat(z);

jacfun = @jacapply;

    function Jv = jacapply(v)
        Jv = cell(Nd,1);
        persistent L U
        
        % Explicit reset proved not to be necessary. Each time this function
        % is redefined for new data, the L and U vars are reset.
        %%if ischar(v), L = []; return, end
        
        if isempty(L)
            L = cell(Nd,1);  U = cell(Nd,1);
            for d = 1:Nd
                [L{d},U{d}] = lu(corjac{d}(z(data(d).idx)));
            end
        end
        
        for d = 1:Nd
            idx = data(d).idx;
            r = zeros(n(d),1);
            r(1) = data(d).Bl*v;
            r(end) = data(d).Br*v;
            Jv{d} = (U{d}\(L{d}\r)) - v(idx);
        end
        Jv = cell2mat(Jv);
    end
end

function [resfun,jacfun] = DAEfuns(u,un,data,dt,lbc,rbc)

idx = data.idx;  M = data.M;
Bl = data.Bl;  Br = data.Br;
uloc = u(idx);
ODEfun = data.resfun;  ODEjac = data.jacfun;

resfun = @fun;
jacfun = @jac;


    function r = fun(z)
        r = ODEfun(z+uloc);
        
        % the true BCs from the formulation of the DAE
        if ~isempty(lbc)
            r(1) = z(1) + uloc(1) - lbc;
        end
        
        if ~isempty(rbc)
            r(end) = z(end) + uloc(end) - rbc;
        end
        
        % without interface interactions:
        r = M*(z+uloc-un(idx)) - dt*r;
        
        % interface conditions        
        if isempty(lbc)
            r(1) = z(1) + uloc(1) - Bl*u;
        end
        
        if isempty(rbc)
            r(end) = z(end) + uloc(end) - Br*u;
        end
        
    end

    function J = jac(z)
        J = ODEjac(z+uloc);
        if ~isempty(lbc)
            J(1,:) = 0;  J(1,1) = 1;
        end
        if ~isempty(rbc)
            J(end,:) = 0;  J(end,end) = 1;
        end
        J = M - dt*J;
        if isempty(lbc)
            J(1,:) = 0;  J(1,1) = 1;
        end
        if isempty(rbc)
            J(end,:) = 0;  J(end,end) = 1;
        end

    end

end


