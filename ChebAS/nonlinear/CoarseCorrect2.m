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
function [ r_er, g_mat,g_p_mat,J_er] = CoarseCorrect2( PUApprox, v,evalf,jacf,tol,tol_c,j)


num_sols = length(v)/length(PUApprox);

v_hat = [];
for i=1:num_sols
    v_hat = [v_hat;PUApprox.Fine2Coarse(v((1:length(PUApprox)) + (i-1)*length(PUApprox)),j)];
end

[g,lg,ug,pg] = ParPreconditionedNewtonForward(v,PUApprox,evalf,jacf,tol);

g_mat.l = lg;
g_mat.u = ug;
g_mat.p = pg;

g_hat = [];

for i=1:num_sols
    g_hat = [g_hat;PUApprox.Fine2Coarse(g((1:length(PUApprox)) + (i-1)*length(PUApprox)),j)];
end

PUApprox.Coarsen();

[g_p,lg_p,ug_p,pg_p] = ParPreconditionedNewtonForward(v_hat,PUApprox,evalf,jacf,tol);

g_p_mat.l = lg_p;
g_p_mat.u = ug_p;
g_p_mat.p = pg_p;


R = g_p-g_hat;
RES = @(u)Residual(u,R,PUApprox,evalf);
JAC = @(u)CoarseASJac(PUApprox,jacf,u,R);

%params = [20,-1,.5,0];
%tol = [1e-5 1e-4];

%[u,~,~,~,~] = nsoldAS(v_hat,RES,JAC,tol,params);

options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',30,'FunctionTolerance',tol_c,'Display','off');
u = fsolve(@(er)sol_and_jac(@(er)RES(er),@(er)JAC(er),er),zeros(size(v_hat)),options);
u = u(:,end);

J_er = JAC(u);

er = u-v_hat;

r_er = [];

for i=1:num_sols
    r_er = [r_er;PUApprox.Coarse2Fine(er((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

PUApprox.Refine();

end

function F = Residual(u,R,PUApprox,evalf)
    F = ParResidualFun(u+R,PUApprox,evalf);
    F = ParResidualInterface(F,u,PUApprox);
end
