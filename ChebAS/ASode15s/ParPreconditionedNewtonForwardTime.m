% INPUT:
%      PUApprox: PUApprox approximation        
%         sol: given solution
%       evalF: residual function which returns Jacobian
%    num_sols: number of solutions
%
% OUTPUT:
%          z: correction of solution
%          J: cell array of local Jacobians
function [z,J] = ParPreconditionedNewtonForwardTime(t,sol,rhs,PUApprox,evalF,num_sols,hinvGak,Jac,M)

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
        rhs_tmp = rhs(loc_sol_step+step(k)+(1:prod(degs)));
        
        sol_loc{k} = [sol_loc{k};tmp];
        rhs_loc{k} = [rhs_loc{k};rhs_tmp];
        
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
    
    [z{k},J{k}] = local_inverse(PUApprox.leafArray{k},sol_loc{k},t,rhs_loc{k},in_border{k},diff{k},evalF,hinvGak,num_sols,length(PUApprox.leafArray{k}),Jac,M{k});
    
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
function [c,J] = local_inverse(approx,sol_k,t,rhs_k,border_k,diff_k,evalF,hinvGak,num_sols,sol_length,Jac,M)

    %The residul is F(sol_k+z_k) 
    %            sol_k(border_k)+z_k(border_k)-B_k*u 
    %            (iterfacing at the zone interface)
    function [F,J] = residual(z)
        
              
        [F] = hinvGak*evalF(approx,t,z+sol_k)-rhs_k-M*z;
        
        z = z+sol_k;  
        
        sol_step = 0;
        
        diff_len = length(diff_k)/num_sols;
        diff_step = 0;
        
        %For each solution, we account for the interfacing.
        %
        % NOTE: I assume each solution has the same length for simplicity
        for i=1:num_sols
            F_s{i} = F(sol_step+(1:sol_length));
            z_loc = z(sol_step+(1:sol_length));
            
            if any(border_k)
                F_s{i}(border_k) = z_loc(border_k) - diff_k(diff_step+(1:diff_len));
            end
            sol_step = sol_step+sol_length;
            diff_step = diff_step + diff_len;
        end
        
        F = cell2mat(F_s');
        
        J = hinvGak*Jac(t,z,approx)-M;
        
        E = eye(sol_length);
        
        for i=1:num_sols
            ind = false(sol_length*num_sols,1);
            ind((i-1)*sol_length+(1:sol_length)) = border_k;
            
            J(ind,:) = zeros(sum(ind),num_sols*sol_length);
            J(ind,(i-1)*sol_length+(1:sol_length)) = E(border_k,:);
        end
        
    end

    options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',1000,'FunctionTolerance',1e-6);

    [s,~,~,~,J] = fsolve(@residual,zeros(size(sol_k)),options);
    
    %opt = [20,-1,.5,0];
    %tol = [1e-2 1e-2];
    %[s, ~, ~, ~,mat_struct] = nsold(zeros(size(sol_k)),@residual,tol,opt);
    
    %L = mat_struct.L;
    %U = mat_struct.U;
   % p = mat_struct.p;
    c = s(:,end);
end


