function [ r_er ] = CoarseCorrect( PUApprox, v,evalf,jacf)


num_sols = length(v)/length(PUApprox);

v_hat = [];
for i=1:num_sols
    v_hat = [v_hat;PUApprox.Fine2Coarse(v((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

r = ParResidual(v,PUApprox,evalf);

r_hat = [];
for i=1:num_sols
    r_hat = [r_hat;PUApprox.Fine2Coarse(r((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

PUApprox.Coarsen();

r_hat = r_hat - ParResidualFun(v_hat,PUApprox,evalf);

params = [20,-1,.5,0];
tol = [1e-4 1e-4];

er = nsoldAS(rand(size(v_hat)),@(er)Residual(er,v_hat,r_hat,PUApprox,evalf),@(er)CoarseASJac(PUApprox,jacf,er,v_hat),tol,params);

r_er = [];

for i=1:num_sols
    r_er = [r_er;PUApprox.Coarse2Fine(er((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

PUApprox.Refine();

end

function F = Residual(er,v,r,PUApprox,evalf)
    F = ParResidualFun(v+er,PUApprox,evalf);
    F = ParResidualInterface(F,er,PUApprox);
    F = F + r;
end

