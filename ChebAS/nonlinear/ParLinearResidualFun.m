% INPUT:     
%      sol: given solution at patches
%      PUApprox: PUApprox approximation   
%      J: Jacobian function for leaf
%
% OUTPUT:
%          z: residual of solution, identity at inner boundary of patches
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [z] = ParLinearResidualFun(sol,PUApprox,J)

%PUApprox.sample(sol);

num_sols = length(sol)/length(PUApprox);

step = zeros(length(PUApprox.leafArray),1);

sol = reshape(sol,length(PUApprox),num_sols);

%Figure out starting index for each patch
for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end

for k=1:length(PUApprox.leafArray)
    
    degs = PUApprox.leafArray{k}.degs;
    
    sol_loc{k} = sol(step(k)+(1:prod(degs)),:);
    
end

%parallel step
for k=1:length(PUApprox.leafArray)
    
    %Assume z is of the form [u1 u2 ... un]
    [z{k}] = local_residual(sol_loc{k},J{k});
    
    z{k} = reshape(z{k},length(PUApprox.leafArray{k}),num_sols);
end

z = cell2mat(z');
z = z(:);

end

% INPUT:
%      approx: leaf approximation        
%       sol_k: given solution
%    border_k: border index for interface
%      diff_k: precomputed interface zone interpolation
%       evalF: residual function which returns Jacobian
%    num_sols: number of solutions
%  sol_length: length of solution
%
% OUTPUT
%           c: correction of solution
%          Jk: local Jocabian
function F = local_residual(sol_k,J)
     
        F = J*sol_k(:);
               
end

