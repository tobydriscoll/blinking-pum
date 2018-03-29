function [ Mat ] = CoarseASMat( Tree,L,B )

if Tree.is_leaf
    LEAVES = {Tree};
else
    LEAVES = Tree.collectLeaves();    
end

Tree.Coarsen();

step = 0;

ii = [];
jj = [];
zz = [];

for k=1:length(LEAVES)
    cdim = LEAVES{k}.cdegs;

    [out_border_c,~,in_border] = FindBoundaryIndex2DSides(cdim,LEAVES{k}.domain,LEAVES{k}.outerbox); 
    
    pointsl = LEAVES{k}.points();
    
    index_n = (1:length(LEAVES{k}))';
    index_n = index_n(in_border);
    
    [iib,jjb,zzb] = Tree.interpMatrixZone_vecs(pointsl(in_border,:));
    
    ii = [ii;index_n(iib)+step];
    jj = [jj;jjb];
    zz = [zz;-zzb];    

    Dx = kron(eye(cdim(2)),diffmat(cdim(1),1,LEAVES{k}.domain(1,:)));
    Dy = kron(diffmat(cdim(2),1,LEAVES{k}.domain(2,:)),eye(cdim(1)));
    
    Dxx = kron(eye(cdim(2)),diffmat(cdim(1),2,LEAVES{k}.domain(1,:)));
    Dyy = kron(diffmat(cdim(2),2,LEAVES{k}.domain(2,:)),eye(cdim(1)));
    
    E = eye(prod(cdim));
    
    OP = L(E,pointsl(:,1),pointsl(:,2),Dx,Dy,Dxx,Dyy);
    
    for i=1:4
        if any(out_border_c{i}) && ~isempty(B{i})
        OP(out_border_c{i},:) = ...
            B{i}(E(out_border_c{i},:),pointsl(out_border_c{i},1),pointsl(out_border_c{i},2),Dx(out_border_c{i},:),Dy(out_border_c{i},:),Dxx(out_border_c{i},:),Dyy(out_border_c{i},:));
        end
    end
    
    OP(in_border,:) = E(in_border,:);
    
    [iid,jjd,zzd] = find(OP);
    
    ii = [ii;iid+step];
    jj = [jj;jjd+step];
    zz = [zz;zzd];
    
    step = step+prod(cdim);
end

Mat = sparse(ii,jj,zz,length(Tree),length(Tree));

Tree.Refine();

end

