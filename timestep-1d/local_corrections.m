function [z,jacfun] = local_corrections(u,un,dt,BC,domain)

% Given data u for the solutions on all domains, find update given by the
% solutionsto local problems using interpolations from u on interfaces.
% Also returns a function that applies the Jacobian of this map at u to any
% vector.

Nd = length(domain);
n = cat(1,domain.n);

z = cell(Nd,1);
opt = [20,-1,.5,0];

daeres = cell(1,Nd);  daejac = cell(1,Nd);
for i = 1:Nd
    if i==1, lbc = BC(1); else lbc = []; end
    if i==Nd, rbc = BC(2); else rbc = []; end
    [daeres{i},daejac{i}] = DAEfuns(u,un,domain(i),dt,lbc,rbc);
    f = nsoldbridge(daeres{i},daejac{i});
    [z{i},ithist] = nsold(zeros(n(i),1),f,[1e-10 1e-10],opt);
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
                [L{d},U{d}] = lu(daejac{d}(z(domain(d).idx)));
            end
        end
        
        Jv = cell(Nd,1);
        for d = 1:Nd
            idx = domain(d).idx;
            r = zeros(n(d),1);
            r(1) = domain(d).Bl*v;
            r(end) = domain(d).Br*v;
            Jv{d} = (U{d}\(L{d}\r)) - v(idx);
        end
        Jv = cell2mat(Jv);
    end
end

function [resfun,jacfun] = DAEfuns(u,un,domain,dt,lbc,rbc)

idx = domain.idx;  M = domain.M;
Bl = domain.Bl;  Br = domain.Br;
uloc = u(idx);
ODEfun = domain.resfun;  ODEjac = domain.jacfun;

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


