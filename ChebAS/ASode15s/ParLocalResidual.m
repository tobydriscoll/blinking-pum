% INPUT:
%      PUApprox: PUApprox approximation        
%         sol: given solution
%       evalF: residual function which returns Jacobian
%    num_sols: number of solutions
%
% OUTPUT:
%          z: correction of solution
%          J: cell array of local Jacobians
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [z] = ParLocalResidual(t,sol,PUApprox,evalF,num_sols)

%PUApprox.sample(sol);

step = zeros(length(PUApprox.leafArray),1);

sol = reshape(sol,length(PUApprox),num_sols);

%Figure out starting index for each patch
for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end

for k=1:length(PUApprox.leafArray)
    
    degs = PUApprox.leafArray{k}.degs;
    
    %This function returns the logical indicies of the gamma and outer
    %boundry interface. Out put is given for all indicies, as well as the
    %indicies along each of the sides
    [~,~,in_border{k},~] = FindBoundaryIndex2DSides(degs,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox);
    
    sol_loc{k} = sol(step(k)+(1:prod(degs)),:);
    
    diff{k} = PUApprox.leafArray{k}.Binterp*sol;
    
end

%parallel step
for k=1:length(PUApprox.leafArray)
    
    [z{k}] = local_residual(PUApprox.leafArray{k},sol_loc{k},t,in_border{k},diff{k},evalF,num_sols);
    
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
function F = local_residual(approx,sol_k,t,border_k,diff_k,evalF,num_sols)
     
        F = evalF(approx,t,sol_k(:));  
        
        sol_length = length(approx);

        F = reshape(F,sol_length,num_sols);
        
        F(border_k,:) = sol_k(border_k,:) - diff_k;
        
        F = F(:);
              
end