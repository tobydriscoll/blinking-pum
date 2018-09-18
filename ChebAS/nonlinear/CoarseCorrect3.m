% INPUT:
%      PUApprox: PUApprox approximation
%             v: given solution
%         evalF: residual function for leaf
%           Jac: jacobian function for leaf
%
% OUTPUT:
%          r_er: correction of solution
%          J_v_pls_er: jacobian of J_hat(v+er)
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [ r_er,J_v_pls_er ] = CoarseCorrect3( PUApprox,Leaf,G,v,evalf,jacf,tol_c)


num_sols = length(v)/length(PUApprox);

v_hat = [];

Len = length(PUApprox);

for i=1:num_sols
    PUApprox.sample( v((1:Len) + (i-1)*Len ) );
    F = PUApprox.evalfGrid(G);
    v_hat = [v_hat;F(:)];
end
Fc = PUchebfun(Leaf);
r = ParResidual(v,PUApprox,evalf);

r_hat = [];

for i=1:num_sols
    PUApprox.sample( r((1:Len) + (i-1)*Len ) );
    F = PUApprox.evalfGrid(G);
    r_hat = [r_hat;F(:)];
end

g = evalf(v_hat,Leaf);
w = g-r_hat; 

params = [20,-1,.5,0];
tol = [1e-5 1e-4];

RES = @(er)Residual(er,v_hat,w,Leaf,evalf);
JAC = @(er)jacf(v_hat+er,Leaf);

%[er,~,~,~,~] = nsoldAS(zeros(size(v_hat)),RES,JAC,[tol_c 10*tol_c],params);

options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',20,'FunctionTolerance',tol_c,'Display','iter');
er = fsolve(@(er)sol_and_jac(@(er)RES(er),@(er)JAC(er),er),zeros(size(v_hat)),options);
er = er(:,end);

J_v_pls_er = JAC(er);

r_er = [];

Len = length(Leaf);

for i=1:num_sols
        
    er_i = er((1:Len) + (i-1)*Len );
    Leaf.sample(er_i);
    
    PUApprox.ChebRoot.sample(@(x)Leaf.evalfGrid(x),true);
    
    F = PUApprox.Getvalues();
    
    r_er = [r_er;F(:)];
    
end

end

function F = Residual(er,v,w,Leaf,evalf)    
    F = evalf(v+er,Leaf) - w;
end

