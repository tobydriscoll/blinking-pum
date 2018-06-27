function [Pz,jacfun] = coarse_correction(u,un,dt,BC,domain,coarse)

% Given data u for the solutions on all domains, find update given by the
% solution on coarse grid, as indicated by the FAS. 
% Also returns a function that applies the Jacobian of this map at u to any
% vector.

Nd = length(domain);
nc = coarse.n;

opt = [20,-1,.5,0];

Ru = coarse.restrict(u);
Run = coarse.restrict(un);
Rr = Ru;

r = cell(Nd,1);
Jr = cell(Nd,1);
for i = 1:Nd
    if i==1, lbc = BC(1); else lbc = []; end
    if i==Nd, rbc = BC(2); else rbc = []; end
    [r{i},Jr{i}] = DAElocal(u,un,domain(i),dt,lbc,rbc);
end
Rr = coarse.restrict(cell2mat(r));

[FRu,J0] = DAE(Ru,Run,dt,BC(1),BC(2),coarse);

    function [y,J] = res(z)
        [y,J] = DAE(z+Ru,Run,dt,BC(1),BC(2),coarse);
        y = y - (FRu-Rr);
    end

[z,ithist,flag] = nsold(zeros(nc,1),@res,[1e-11 1e-11],opt);
if flag > 0, warning('nsold did not finish well in coarse correction'), end

Pz = coarse.prolong(z); 

jacapply('clear')
jacfun = @jacapply;

    function Jv = jacapply(v)
        persistent L U
        
        % Not clear whether explicit reset is necessary. 
        if ischar(v), L = []; return, end
        
        if isempty(L)
            [~,J] = DAE(z+Ru,Run,dt,BC(1),BC(2),coarse);
            [L,U] = lu(J);
        end
        
        Rv = coarse.restrict(v);
        drdu_v = zeros(sum(cat(1,domain.n)),1);
        for d = 1:Nd
            drdu_v(domain(d).idx) = Jr{d}*v(domain(d).idx);
        end
        
        
        w = -Rv + U\(L\(J0*Rv-coarse.restrict(drdu_v)));
        Jv = coarse.prolong(w);
        
     end
end

function [r,J] = DAElocal(u,un,data,dt,lbc,rbc)

idx = data.idx;  M = data.M;
ODEfun = data.resfun;  ODEjac = data.jacfun;

ui = u(idx);

r = ODEfun(ui);

% the true BCs from the formulation of the DAE
if ~isempty(lbc)
    r(1) = ui(1) - lbc;
end
       
if ~isempty(rbc)
    r(end) = ui(end) - rbc;
end

% without interface interactions:
r = M*(ui-un(idx)) - dt*r;
               
% interface conditions
if isempty(lbc)
    r(1) = ui(1) - data.Bl*u;
end
        
if isempty(rbc)
    r(end) = ui(end) - data.Br*u;
end
        
J = ODEjac(ui);
        
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


function [r,J] = DAE(u,un,dt,lbc,rbc,data)

M = data.M;
ODEfun = data.resfun;  ODEjac = data.jacfun;

r = ODEfun(u);

r(1) = u(1) - lbc;
r(end) = u(end) - rbc;

r = M*(u-un) - dt*r;
        
J = ODEjac(u);
J(1,:) = 0;  J(1,1) = 1;
J(end,:) = 0;  J(end,end) = 1;

J = M - dt*J;

end



