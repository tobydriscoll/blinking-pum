% INPUT:
%             t: time in the ODE
%           sol: given solution
%           rhs: rhs term. Typically from time integration method
% PUApproxArray: cell array of Trees for each solution. Each tree has
%                interface indicies for each leaf, as well as boundary
%                info if solution is to be 'packed'.
%
%                It is expected that the patch structure is the same
%                between different solutions. Patches with respective
%                solutions can have different number of unknowns though.
%     NonLinOps: cell array of Nonlin ops for each leaf (for our problem,
%                it is the blink objects).
%       hinvGak: this is the time step term used within ode15s
%             M: cell array of mass matrix
%
% OUTPUT:
%          z: correction of solution
%          l,u,p: cell arrays of local Jacobian factorizations
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [z,l,u,p,resnorm] = SNKresidual(t,sol,rhs,PUApproxArray,NonLinOps,hinvGak,M,alpha,use_par)

rhs_zero = rhs==0;

if rhs_zero
    rhs = zeros(size(sol));
end

num_sols = length(PUApproxArray);
num_leaves = length(PUApproxArray{1}.leafArray);
sol_lengths = zeros(num_sols,1);
for i=1:num_sols
    sol_lengths(i) = length(PUApproxArray{i});
end

if ~iscell(M)
    M = repmat({0},1,num_leaves);
end

sol = mat2cell(sol,sol_lengths);
sol_unpacked = sol;

for i=1:num_sols
    
    if PUApproxArray{i}.is_packed
        %pull the boundry info for the packed functions. The boundary
        %info is stored within the tree itself.
        sol_unpacked{i} = PUApproxArray{i}.Getunpackedvalues(sol{i});
    else
        sol_unpacked{i} = sol{i};
    end
end

[ sol_loc,lens ] = unpackPUvecs(cell2mat(sol),PUApproxArray);
rhs_loc = unpackPUvecs(rhs,PUApproxArray);
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
if use_par
	parfor k=1:num_leaves
        [z_loc{k},l{k},u{k},p{k},resnorm{k}] = local_inverse(sol_loc{k},t,rhs_loc{k},diff{k},border{k},NonLinOps{k},hinvGak,M{k},lens{k},alpha);
	end
else
	for k=1:num_leaves
        [z_loc{k},l{k},u{k},p{k},resnorm{k}] = local_inverse(sol_loc{k},t,rhs_loc{k},diff{k},border{k},NonLinOps{k},hinvGak,M{k},lens{k},alpha);
	end
end

z = packPUvecs(z_loc,PUApproxArray);
resnorm = norm([resnorm{:}]);

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
function [c,l,u,p,resnorm] = local_inverse(sol_k,t,rhs_k,diff_k,border_k,NonLinOps_k,hinvGak,M,lens_k,alpha)

%The residul is F(sol_k+z_k)
%            sol_k(border_k)+z_k(border_k)-B_k*u
%            (iterfacing at the zone interface)
    function F = residual(z)
        num_sols = length(lens_k);
        F = hinvGak*NonLinOps_k.timederiv(t,z+sol_k)-M*(z+rhs_k);
        F = mat2cell(F,lens_k);
        z = mat2cell(z,lens_k);
        sol_k_c =  mat2cell(sol_k,lens_k);
        for i=1:num_sols
            F{i}(border_k{i}) = alpha*(z{i}(border_k{i}) +sol_k_c{i}(border_k{i}) - diff_k{i});
        end
        F = cell2mat(F);
    end

    function J = jac_fun(z)
        num_sols = length(lens_k);
        J = hinvGak*NonLinOps_k.jac(t,z+sol_k)-M;
        index = 0;
        
        %This is supposed to account for the interfacing
        for i=1:num_sols
            E = eye(lens_k(i));
            total_length = sum(lens_k);
            ind = false(total_length,1);
            local_ind = index+(1:lens_k(i));
            ind(local_ind) = border_k{i};
            
            J(ind,:) = zeros(sum(ind),total_length);
            J(ind,local_ind) = alpha*E(border_k{i},:);
            index = index+lens_k(i);
        end
    end

% options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',15,'Display','off');
% [c,~,~,~,~] = fsolve(@(u)sol_and_jac(@residual,@jac_fun,u),zeros(numel(sol_k),1),options);
% c = c(:,end);
% J = jac_fun(c);

%[l,u,p] = lu(J,'vector');

params = [30,-1,.5,0];
tol = [1e-6 1e-6];

[c,l,u,p,jac,ierr,stat] = nsoldAS(zeros(numel(sol_k),1),@residual,@jac_fun,tol,params);
resnorm = stat(end,2);

%fprintf('norm(c) = %.3e\n',norm(c))
% a=1;
%c = s(:,end);
%if norm(c)==0, keyboard, end
end


