% INPUT:     
%      sol: given solution at patches
%      sol2: second given solution 
%      PUApprox: PUApprox approximation   
%      jac_f(x,leaf): function that returns jacobian given local solution x and
%      patch leaf.
%
% OUTPUT:
%          Mat: sparse matrix used for Jacobian of Coarse correction.
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [ Mat ] = CoarseInterfaceInterp(PUApprox,num_sols)
%assume sol is the correct coarse length

PUApprox.Coarsen();

step = zeros(length(PUApprox.leafArray),1);

%Figure out starting index for each patch
for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + length(PUApprox.leafArray{k-1});
end

ii = [];
jj = [];
zz = [];

for k=1:length(PUApprox.leafArray)
    
    degs = PUApprox.leafArray{k}.cdegs;
    
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
    
end

Mat = sparse(ii,jj,zz,num_sols*length(PUApprox),num_sols*length(PUApprox));

PUApprox.Refine();

end



