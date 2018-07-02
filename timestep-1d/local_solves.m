function [z,jacfun] = local_solves(u,un,dt,BC,domain,data)

% Given data u for the solutions on all domains, find the solutions given by
% local solves using interpolations from u on interfaces.
% Also returns a function that applies the Jacobian of this map at u to any
% vector.

Nd = length(domain);
n = cat(1,domain.n);

z = cell(Nd,1);
opt = [20,-1,.5,0];

DAEres = cell(1,Nd);  DAEjac = cell(1,Nd);
for i = 1:Nd
    if i==1, lbc = BC(1); else lbc = []; end
    if i==Nd, rbc = BC(2); else rbc = []; end
    [DAEres{i},DAEjac{i}] = DAEfuns(u,un,data(i),dt,lbc,rbc);
    f = nsoldbridge(DAEres{i},DAEjac{i});
    [z{i},ithist] = nsold(u(data(i).idx),f,[1e-11 1e-11],opt);
end
z = cell2mat(z);

jacfun = @jacapply;

    function Jv = jacapply(v)
        persistent L U
        
        % Not clear whether explicit reset is necessary. 
        if ischar(v), L = []; return, end
        
        if isempty(L)
            L = cell(Nd,1);  U = cell(Nd,1);
            for d = 1:Nd
                [L{d},U{d}] = lu(DAEjac{d}(z(data(d).idx)));
            end
        end
        
        Jv = cell(Nd,1);
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
ODEfun = data.resfun;  ODEjac = data.jacfun;

resfun = @fun;
jacfun = @jac;


    function r = fun(z)
        r = ODEfun(z);
        
        % the true BCs from the formulation of the DAE
        if ~isempty(lbc)
            r(1) = z(1) - lbc;
        end
        
        if ~isempty(rbc)
            r(end) = z(end) - rbc;
        end
 
        % without interface interactions:
        r = M*(z-un(idx)) - dt*r;
               
        % interface conditions
        if isempty(lbc)
            r(1) = z(1) - Bl*u;
        end
        
        if isempty(rbc)
            r(end) = z(end) - Br*u;
        end
        
    end

    function J = jac(z)
        J = ODEjac(z);
        
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


