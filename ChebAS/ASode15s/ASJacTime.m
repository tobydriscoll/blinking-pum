% CoarseASJac
% This method computes the Jacobian used in the nonlinear solve of the
% coarse correction in the SNK two level method.
%
% INPUT:     
%      PUApprox: PUApprox approximation   
%
%      Jac: Jacobian function Jac(x,p) for solution x and local
%             approximation p. Jac(x,p) evaluates the residual on the domain
%             of p.
%
%       sol: given solution at patches
%
%       sol2: second given solution (this can probably be removed)
%
% OUTPUT:
%          Mat: sparse matrix used for Jacobian of coarse correction.
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [Mat,J_sol] = ASJacTime(PUApproxArray,NonlinOp,M,dt,t,sol,alpha,rhs)
%assume sol is the correct coarse length

J_sol = [];

find_J_rhs = false;

if nargin>7
    find_J_rhs = true;
    
    rhs_loc = unpackPUvecs(rhs,PUApproxArray);
end
    
if ~iscell(PUApproxArray)
    PUApproxArray = {PUApproxArray};
end

num_sols = length(PUApproxArray);

num_leaves = length(PUApproxArray{1}.leafArray);

sol_lengths = zeros(num_sols,1);

for i=1:num_sols
    sol_lengths(i) = length(PUApproxArray{i});
end

ii = [];
jj = [];
zz = [];

[sol_loc,len_loc] = unpackPUvecs(sol,PUApproxArray);

sol_index = 0;

sub_ind = cell(1,num_sols);



%We need to ignore outer boundary data on patches that are packed, i.e.
%solutions with assumed dirichlet data
for i=1:num_sols
    
    sub_ind{i} = true(0,1);
    
    for k=1:num_leaves
    
        loc_leaf = PUApproxArray{i}.leafArray{k};
        
        if PUApproxArray{i}.is_packed
            sub_ind{i} = [sub_ind{i};~loc_leaf.outer_boundary];
        else
            sub_ind{i} = [sub_ind{i};true(length(loc_leaf),1)];
        end
    end
end

for i=1:num_sols
    %    Figure out sparse indicies of interpolation matrix
    leaf_index = 0;
    
    for k=1:num_leaves
        
        loc_leaf = PUApproxArray{i}.leafArray{k};
        
        B = loc_leaf.Binterp;
        
        B = B(:,sub_ind{i});
        
        [iib,jjb,zzb] = find(B);
        

        
        index_n = (1:length(loc_leaf))';
        index_n = index_n(loc_leaf.inner_boundary);
        
        %order is solution then patch
        
        %first index is to the inner boundary of the k_th patch
        %in the i_th solution.
        iib = index_n(iib)+sol_index+leaf_index;
        
        %second index is from the j_th solution.
        jjb = jjb + sol_index;
        
        %These are just the values of the matrix
        zzb = -zzb;
        
        ii = [ii;iib];
        jj = [jj;jjb];
        zz = [zz;zzb];
        
        leaf_index = leaf_index + length(loc_leaf);
        
        border{k}{i} = loc_leaf.inner_boundary;
    end
    
    sol_index = sol_index+sol_lengths(i);
    
end


loc_leaf_index = zeros(num_leaves,num_sols);

%For [u_1 u_2;v_1 v_2] determine the starting index for each local
%solution i.e. loc_leaf_index(2,2) returns the starting index of v_2
%within [v_1 v_2]
for k=2:num_leaves
   loc_leaf_index(k,:) = loc_leaf_index(k-1,:)+len_loc{k-1};
end


for k=1:num_leaves
    
    J = loc_jacobian(NonlinOp{k},sol_loc{k},t,dt,M{k},len_loc{k},border{k});
    
    if find_J_rhs
        J_sol{k} = J*rhs_loc{k};
    end
    
    [iid,jjd,zzd] = find(J);
    
    shift_idd = iid;
    shift_jdd = jjd;
    
    loc_sol_index = 0;
    global_sol_index = 0;
    
    for i=1:num_sols
         
        loc_sol_length = len_loc{k}(i);
        
        %find the indicies of v_i
        ind = (iid >= (loc_sol_index + 1)) & (iid <= (loc_sol_index + loc_sol_length)); 
        
        %take the indicies of v_i to 1:length(v_i)
        shift_idd(ind) = iid(ind) - loc_sol_index;
        
        %adding global_sol_index places the v_i withen v in [u;v]
        % adding loc_leaf_index(k,i) places v_i in the i_th patch of
        % [v_1;v_2;...v_n].
        shift_idd(ind) = shift_idd(ind)+global_sol_index+loc_leaf_index(k,i);
        
        
        %do the same for j.
        ind = (jjd >= (loc_sol_index + 1)) & (jjd <= (loc_sol_index + loc_sol_length)); 
        
        shift_jdd(ind) = jjd(ind) - loc_sol_index;
        
        shift_jdd(ind) = shift_jdd(ind)+global_sol_index+loc_leaf_index(k,i);
        
        
        %update local solution and global solution starting index
        loc_sol_index = loc_sol_index+loc_sol_length;
        
        global_sol_index = global_sol_index+sol_lengths(i);
        
    end
    
    ii = [ii;shift_idd];
    jj = [jj;shift_jdd];
    zz = [zz;zzd];
    
end    

Mat = sparse(ii,jj,zz,sum(sol_lengths),sum(sol_lengths));

if find_J_rhs
    J_sol = packPUvecs(J_sol,PUApproxArray);
end

end

function J = loc_jacobian(NonLinOps_k,sol_k,t,dt,M,lens_k,border_k)
        
        num_sols = length(lens_k);
        
        J = dt*NonLinOps_k.jac(t,sol_k)-M;

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



