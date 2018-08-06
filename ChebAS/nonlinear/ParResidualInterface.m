% INPUT:     
%      sol: given solution at patches
%      sol2: second given solution 
%      PUApprox: PUApprox approximation   
%      evalF: residual function which returns Jacobian
%
% OUTPUT:
%          z: returns sol, with sol-sol2 of inner boundary of patches
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [z] = ParResidualInterface(sol,sol2,PUApprox)

%PUApprox.sample(sol);

num_sols = length(sol)/length(PUApprox);

step = zeros(length(PUApprox.leafArray),1);

sol = reshape(sol,length(PUApprox),num_sols);
sol2 = reshape(sol2,length(PUApprox),num_sols);

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
    
    if ~PUApprox.iscoarse
        diff{k} = PUApprox.leafArray{k}.Binterp*sol2;
    else
        diff{k} = PUApprox.leafArray{k}.CBinterp*sol2;
    end
    
end

%parallel step
for k=1:length(PUApprox.leafArray)
    
    %Assume z is of the form [u1 u2 ... un]
    [z{k}] = local_residual(sol_loc{k},in_border{k},diff{k});
    
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
function F = local_residual(sol_k,border_k,diff_k)
        
        F = sol_k;
        
        F(border_k,:) = sol_k(border_k,:) - diff_k;
              
end

