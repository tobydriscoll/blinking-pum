% INPUT:
%      PUApprox: PUApprox approximation        
%         sol: given solution
%       evalF: residual function which returns Jacobian
%    num_sols: number of solutions
%
% OUTPUT:
%          z: correction of solution
%          J: cell array of local Jacobians
function [z] = ParLocalResidual(t,sol,PUApprox,evalF,num_sols)

%PUApprox.sample(sol);

step = zeros(length(PUApprox.leafArray),1);

%Figure out starting index for each patch
for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + num_sols*length(PUApprox.leafArray{k-1});
end

for k=1:length(PUApprox.leafArray)
    
    degs = PUApprox.leafArray{k}.degs;
    
    sol_loc{k} = [];
    rhs_loc{k} = [];
    diff{k} = [];
    
    %This function returns the logical indicies of the gamma and outer
    %boundry interface. Out put is given for all indicies, as well as the
    %indicies along each of the sides
    [~,~,in_border{k},~] = FindBoundaryIndex2DSides(degs,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox);
    
    loc_sol_step = 0;
    sol_step = 0;
    
    %Some care is given when the PDE has multiple functions to solve for.
    %Here interfacing is precomputed for each solution.
    %
    % NOTE: I assume each solution has the same length for simplicity
    for j=1:num_sols
        tmp = sol(loc_sol_step+step(k)+(1:prod(degs)));
        
        sol_loc{k} = [sol_loc{k};tmp];
        
        %This computes solution interpolated on the zone interface of the
        %patch.
        
        if length(PUApprox.leafArray)>1
            diff{k} = [diff{k};PUApprox.leafArray{k}.Binterp*sol(sol_step+(1:length(PUApprox)))];
            loc_sol_step = loc_sol_step + prod(degs);
            sol_step = sol_step + length(PUApprox);
        end
    end
    
end

%parallel step
for k=1:length(PUApprox.leafArray)
    
    [z{k}] = local_residual(PUApprox.leafArray{k},sol_loc{k},t,in_border{k},diff{k},evalF,num_sols,length(PUApprox.leafArray{k}));
    
end

z = cell2mat(z');

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
function F = local_residual(approx,sol_k,t,border_k,diff_k,evalF,num_sols,sol_length)
     
        F = evalF(approx,t,sol_k);  
        
        sol_step = 0;
        
        diff_len = length(diff_k)/num_sols;
        diff_step = 0;
        
        %For each solution, we account for the interfacing.
        %
        % NOTE: I assume each solution has the same length for simplicity
        for i=1:num_sols
            F_s{i} = F(sol_step+(1:sol_length));
            s_loc = sol_k(sol_step+(1:sol_length));
            
            if any(border_k)
                F_s{i}(border_k) = s_loc(border_k) - diff_k(diff_step+(1:diff_len));
            end
            sol_step = sol_step+sol_length;
            diff_step = diff_step + diff_len;
        end
        
        F = cell2mat(F_s');
        
        
end