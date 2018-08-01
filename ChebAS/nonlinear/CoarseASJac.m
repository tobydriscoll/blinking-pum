% This method sets up the linear operators to be used in the patches for
% the theta method.
%
%   input:
%    Tree: Tree that has the current time step stored
%       L: Operator of PDE u_t = L u
%       B: cell array of boundary condition operator in order NORTH SOUTH EAST WEST
%   theta: parameter used in theta method
%    dt,t: time step and current time
%  output:
%     MAT: returns the coarse matrix for the theta method
function [ Mat ] = CoarseASJac(PUApprox,jac_f,sol,sol2)
%assume sol is the correct coarse length
num_sols = length(sol)/length(PUApprox);

step = zeros(length(PUApprox.leafArray),1);

sol = reshape(sol,length(PUApprox),num_sols);
sol2 = reshape(sol2,length(PUApprox),num_sols);

%Figure out starting index for each patch
for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end

ii = [];
jj = [];
zz = [];

for k=1:length(PUApprox.leafArray)
    
    degs = PUApprox.leafArray{k}.cdegs;
    
    sol_loc = sol(step(k)+(1:prod(degs)),:);
    sol2_loc = sol2(step(k)+(1:prod(degs)),:);
    
    %Figure out indicies of boundary and interface points
    [~,~,in_border,~]  = FindBoundaryIndex2DSides(degs,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox);
        
    index_n = (1:length(PUApprox.leafArray{k}))';
    index_n = index_n(in_border);
    
%    Figure out sparse indicies of interpolation matrix
    [iib,jjb,zzb] = find(PUApprox.leafArray{k}.CBinterp);
    % add them to the matrix
    % Note that ibb is the indicies of the interface border points.
    % This is why we use index_n(iib) here.
    
    for i=1:num_sols
        %current solution, then patch
        ii = [ii;index_n(iib)+(i-1)*length(PUApprox)+step(k)];
        jj = [jj;jjb+(i-1)*length(PUApprox)];
        zz = [zz;-zzb];
    end
     
    J = jac_f(sol_loc(:)+sol2_loc(:),PUApprox.leafArray{k});
    
    
    sol_length = prod(degs);
    
    E = eye(sol_length);
    
    for i=1:num_sols
        ind = false(sol_length*num_sols,1);
        ind((i-1)*sol_length+(1:sol_length)) = in_border;
        
        J(ind,:) = zeros(sum(ind),num_sols*sol_length);
        J(ind,(i-1)*sol_length+(1:sol_length)) = E(in_border,:);
    end
    
    [iid,jjd,zzd] = find(J);

    shift_idd = iid;
    shift_jdd = jjd;

    for i=2:num_sols       
        
        ind = (iid >= ((i-1)*sol_length + 1)) & (iid <= ((i-1)*sol_length + sol_length)); 
        
        shift_idd(ind) = iid(ind) - (i-1)*sol_length;
        shift_idd(ind) = shift_idd(ind) + (i-1)*length(PUApprox);
        
        ind = (jjd >= ((i-1)*sol_length + 1)) & (jjd <= ((i-1)*sol_length + sol_length)); 
        
        shift_jdd(ind) = jjd(ind) - (i-1)*sol_length;
        shift_jdd(ind) = shift_jdd(ind) + (i-1)*length(PUApprox);
        
    end
    
    ii = [ii;shift_idd+step(k)];
    jj = [jj;shift_jdd+step(k)];
    zz = [zz;zzd];
    
end

Mat = sparse(ii,jj,zz,num_sols*length(PUApprox),num_sols*length(PUApprox));

end



