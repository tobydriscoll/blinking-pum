% INPUT:
%       PUApprox: Cell Array of PUApprox approximation
%
%           sol: given solution
%
%           rhs: rhs term. Typically from time integration method
%
% PUApproxArray: cell array of Trees for each solution. Each tree has
%                interface indicies for each leaf, as well as boundary
%                info if solution is to be 'packed'.
%                 
%                It is expected that the patch structure is the same
%                between different solutions. Patches with respective
%                solutions can have different number of unknowns though.
%
%     NonLinOps: cell array of Nonlin ops for each leaf (for our problem,
%                it is the blink objects).
%
% OUTPUT:
%          z: correction of solution
%          J: cell array of local Jacobians
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [z,l,u,p] = SNK_forward_eval(sol,PUApproxArray,NonLinOps,rhs,params)

if ~iscell(PUApproxArray)
    PUApproxArray = {PUApproxArray};
end

if nargin<4
    rhs = 0;
end

if nargin<5
 params = [];
end

num_sols = length(PUApproxArray);

num_leaves = length(PUApproxArray{1}.leafArray);

sol_lengths = zeros(num_sols,1);

for i=1:num_sols
    sol_lengths(i) = length(PUApproxArray{i});
end

sol = mat2cell(sol,sol_lengths);

sol_unpacked = sol;

for i=1:num_sols
    
    if PUApproxArray{i}.is_packed
        
        setBoundary(NonLinOps,PUApproxArray);
        %pull the boundry info for the packed functions. The boundary
        %info is stored within the tree itself.
        sol_unpacked{i} = PUApproxArray{i}.Getunpackedvalues(sol{i});
        
    else
        
        sol_unpacked{i} = sol{i};
    
    end
    
end

%Take [u1;u2;v1;v2] to {[u1;u2],[v1;v2]}
[ sol_loc,lens ] = unpackPUvecs(cell2mat(sol),PUApproxArray);

if rhs~=0
    rhs_loc = unpackPUvecs(rhs,PUApproxArray);
end

start_index = zeros(num_sols,1);

diff = cell(num_leaves,num_sols);
border = cell(num_leaves,num_sols);

for k=1:num_leaves

    for i=1:num_sols
        %This will be (interface length)*num_sols
        diff{k}{i} = PUApproxArray{i}.leafArray{k}.Binterp*sol_unpacked{i};
        border{k}{i} = PUApproxArray{i}.leafArray{k}.inner_boundary;
    end
    
    start_index = start_index + lens{k};
    
end


%parallel step
for k=1:num_leaves
    
    if rhs~=0
        [z_loc{k},l{k},u{k},p{k}] = local_inverse(sol_loc{k},rhs_loc{k},diff{k},border{k},NonLinOps{k},lens{k},params);
    else
        [z_loc{k},l{k},u{k},p{k}] = local_inverse(sol_loc{k},0,diff{k},border{k},NonLinOps{k},lens{k},params);
    end
    
end

%Take {[u1;u2],[v1;v2]} to [u1;u2;v1;v2]
z = packPUvecs(z_loc,PUApproxArray);

end

% INPUT:
%       sol_k: given solution
%           t: current time
%       rhs_k: local right hand side
%      diff_k: precomputed interface zone interpolation
%    border_k: border index for interface for each solution
%   NonLinOps: 
%       evalF: residual function which returns Jacobian
%    num_sols: number of solutions
%  sol_length: length of solution
%
% OUTPUT
%           c: correction of solution
%          Jk: local Jocabian
function [c,l,u,p] = local_inverse(sol_k,rhs_k,diff_k,border_k,NonLinOps_k,lens_k,params)

%The residul is F(sol_k+z_k)
%            sol_k(border_k)+z_k(border_k)-B_k*u
%            (iterfacing at the zone interface)
    function [F] = residual(z)
        
        num_sols = length(lens_k);
        
        F = NonLinOps_k.residual(z+sol_k,params)+rhs_k;
        
        F = mat2cell(F,lens_k);
        
        z = mat2cell(z,lens_k);
        
        sol_k_c =  mat2cell(sol_k,lens_k);
        
        for i=1:num_sols
            F{i}(border_k{i}) = z{i}(border_k{i}) + sol_k_c{i}(border_k{i}) - diff_k{i};
        end
        
       F = cell2mat(F);
        
    end

    function J = jac_fun(z)
        
        num_sols = length(lens_k);
        
        J = NonLinOps_k.jac(z+sol_k,params);

        index = 0;
        
        %This is supposed to account for the interfacing
        for i=1:num_sols
            
            E = eye(lens_k(i));
            
            total_length = sum(lens_k);
            
            ind = false(total_length,1);
            
            local_ind = index+(1:lens_k(i));
            
            ind(local_ind) = border_k{i};
            
            J(ind,:) = zeros(sum(ind),total_length);
            J(ind,local_ind) = E(border_k{i},:);
            
            index = index+lens_k(i);
        end
    end


options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',1000,'FunctionTolerance',1e-4,'Display','iter');
[c,~,~,~,J] = fsolve(@(u)sol_and_jac(@residual,@jac_fun,u),zeros(numel(sol_k),1),options);
c = c(:,end);

[l,u,p] = lu(J,'vector');
end

function setBoundary(NonLinOps,PUApproxArray)
    for i=1:length(PUApproxArray{1}.leafArray)
        NonLinOps{i}.SetBoundaryVals();
    end
end


