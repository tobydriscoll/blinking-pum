% INPUT:
%      PUApprox: PUApprox approximation        
%            w: solution to be applied to Jacobian
%               patches of f(v+e)
%      Jac_hat: f_hat'(v+e)
%          FJv: cell array of f'(v_j)
%       FJv_hat: matrix of f_hat'(v_hat)
%
% NOTE x is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
function [ y ] = LinearGlobalCoarseCorrect( PUApprox, Leaf, G, w,Jac_hat)


num_sols = length(w)/length(PUApprox);

w_hat = [];

Len = length(PUApprox);

for i=1:num_sols
    PUApprox.sample( w((1:Len) + (i-1)*Len ) );
    F = PUApprox.evalfGrid(G);
    w_hat = [w_hat;F(:)];
end

y_hat = Jac_hat.U\Jac_hat.L\w_hat(Jac_hat.p);

Len = length(Leaf);

y = [];

for i=1:num_sols
        
    y_hat_i = y_hat((1:Len) + (i-1)*Len );
    Leaf.sample(y_hat_i);
    
    PUApprox.ChebRoot.sample(@(x)Leaf.evalfGrid(x),true);
    
    F = PUApprox.Getvalues();
    
    y = [y;F(:)];
    
end

end
