% INPUT:
%      PUApprox: Cell Array of PUApprox approximation
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
function [z,l,u,p] = ParPreconditionedNewtonForwardTime(t,sol,rhs,PUApproxArray,NonLinOps,hinvGak,Jac,M)

num_sols = length(PUApproxArray);

sol_lengths = zeros(num_sols,1);

for i=1:num_sols
    sol_lengths(i) = length(PUApproxArray{i});
end

sol = mat2cell(sol,sol_lengths);
rhs = mat2cell(rhs,sol_lengths);

sol_unpacked = sol;

for i=1:num_sols
    
    if PUApproxArray{i}.is_packed
        
        PUApproxArray{i}.Setvalues(sol{i});
        PUApproxArray{i}.unpack();
        sol_unpacked{i} = PUApproxArray{i}.Getvalues();
        PUApproxArray{i}.pack();
        
    else
        
        sol_unpacked{i} = sol{i};
    
    end
    
end

start_index = zeros(num_sols,1);

for k=1:length(PUApprox.leafArray)
    
    sol_loc{k} = [];
    rhs_loc{k} = [];
    lens{k} = [];
    PUApproxArray{i}.leafArray{k};
    
    border{k} = cellarray(1,3);
    
    for i=1:num_sols
        len = length(PUApproxArray{i}.leafArray{k});
        sol_loc{k} = [sol_loc{k};sol{i}(start_index(i):start_index(i)+len)];
        rhs_loc{k} = [rhs_loc{k};rhs{i}(start_index(i):start_index(i)+len)];
        
        lens = [lens len];

        border{k}{i} = PUApproxArray{i}.inner_boundary;
        
    end
    
    for i=1:num_sols
        %This will be (interface length)*num_sols
        diff{k}{i} = PUApproxArray{i}.leafArray{k}.Binterp*sol_unpacked{i};
    end
end

%parallel step

leafs = PUApprox.leafArray;

for k=1:length(leafs)
    
    [z{k},l{k},u{k},p{k}] = local_inverse(leafs{k},sol_loc{k},t,rhs_loc{k},diff{k},border{k},NonLinOps{k},hinvGak,num_sols,Jac,M{k},lens);
    
    z{k} = reshape(z{k},length(leafs{k}),num_sols);
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
function [c,l,u,p] = local_inverse(sol_k,t,rhs_k,diff_k,border_k,NonLinOps_k,hinvGak,num_sols,M,lens_k)

%The residul is F(sol_k+z_k)
%            sol_k(border_k)+z_k(border_k)-B_k*u
%            (iterfacing at the zone interface)
    function [F] = residual(z)
        
        F = hinvGak*NonLinOps_k(t,z+sol_k)+rhs_k(:)-M*z;
        
        F = mat2cell(F,lens_k);
        
        z = mat2cell(F,lens_k);
        
        sol_k =  mat2cell(F,lens_k);
        
        for i=1:num_sols
            F{i}(border_k{i}) = z{i}(border_k{i}) + sol_k{i}(border_k{i}) - diff_k{i};
        end
        
        F = cell2mat(F);
        
    end

    function J = jac_fun(z)
        
        J = NonLinOps_k.jac(t,z+sol_k)-M;

        index = 0;
        
        for i=1:num_sols
            
            E = eye(lens_k(i));
            
            total_length = sum(lens_k);
            
            ind = false(total_length,1);
            
            ind(index+(1:lens_k(i))) = border_k{i};
            
            J(ind,:) = zeros(sum(ind),total_length);
            J(ind,index+(1:lens_k(i))) = E(border_k{i},:);
            
        end
    end


options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',1000,'FunctionTolerance',1e-4,'Display','off');
[c,~,~,~,~] = fsolve(@(u)sol_and_jac(@residual,@jac_fun,u),zeros(numel(sol_k),1),options);
c = c(:,end);

J = jac_fun(c);

[l,u,p] = lu(J,'vector');

%params = [20,-1,.5,0];
%tol = [1e-4 1e-4];
%[c,l,u,p] = nsoldAS(zeros(numel(sol_k),1),@residual,@jac_fun,tol,params);

%c = s(:,end);
end


