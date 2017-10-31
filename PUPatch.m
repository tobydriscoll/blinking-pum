classdef PUPatch<Patch
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    % To do. figure out way to deal with weights. Maybe have seperate
    % object for weights? IDK.
    properties
        children
        splitting_dim
        overlap_in
        index = [];
    end
    
    properties (Constant)
        weights = PUWeights();
    end
    
    methods
        
        function obj = PUPatch(region,zone,overlap_in,cheb_length,children,splitting_dim,index)
            obj.outerbox = children{1}.outerbox;
            obj.region = region;
            obj.zone = zone;
            obj.domain = region;
            [obj.dim,~] = size(obj.domain);
            obj.overlap_in = overlap_in;
            obj.cheb_length = cheb_length;
            obj.children = children;
            obj.splitting_dim = splitting_dim;
            obj.is_leaf = false;
            obj.is_refined = false;
            obj.children{1}.index = [index 1];
            obj.children{2}.index = [index 2];
        end
        
        function ln = length(obj)
            ln = length(obj.children{1})+length(obj.children{2});
        end
        
        function grids = leafGrids(obj)
            gridsl = leafGrids(obj.children{1});
            gridsr = leafGrids(obj.children{2});
            if obj.children{1}.is_leaf && obj.children{2}.is_leaf
                grids = {gridsl gridsr};
            elseif obj.children{1}.is_leaf && ~ obj.children{2}.is_leaf
                grids = {gridsl gridsr{:}};
            elseif obj.children{2}.is_leaf
                grids = {gridsl{:} gridsr};
            else
                grids = {gridsl{:} gridsr{:}};
            end
        end
        
        %Implement this later.
        function pts = points(obj)
            pts = [obj.children{1}.points();obj.children{2}.points()];
        end
        
        function is_refined = PUsplit(obj,Max,set_vals)
            
            if nargin ==2
                set_vals = false;
            end
            
            is_refined = true;
            is_geometric_refined = true;
            
            for k=1:2
                
                if ~obj.children{k}.is_refined
                    if obj.children{k}.is_leaf
                        obj.children{k} = obj.children{k}.splitleaf(Max,set_vals);
                        is_refined = is_refined && obj.children{k}.is_refined;
                        is_geometric_refined = is_geometric_refined && obj.children{k}.is_geometric_refined;
                    else
                        temp = obj.children{k}.PUsplit(Max,set_vals);
                        is_refined = is_refined && temp;
                        is_geometric_refined = is_geometric_refined && obj.children{k}.is_geometric_refined;
                    end
                end
            end
            
            obj.domain = [obj.children{1}.region(:,1) obj.children{2}.region(:,2)];
            
            obj.cheb_length = obj.children{1}.cheb_length+obj.children{2}.cheb_length;
            obj.overlap_in = [obj.children{2}.domain(obj.splitting_dim,1), obj.children{1}.domain(obj.splitting_dim,2)];
            
            obj.is_refined = is_refined;
            obj.is_geometric_refined = is_geometric_refined;
        end
        
        function [Max] = sample(obj,f)
            if ~isnumeric(f)
                Max1 = obj.children{1}.sample(f);
                Max2 = obj.children{2}.sample(f);
            else
                Max1 = obj.children{1}.sample(f(1:length(obj.children{1})));
                Max2 = obj.children{2}.sample(f(length(obj.children{1})+1:end));
            end
            Max = max(Max1,Max2);
        end
        
        
        %  findIndex(obj,index)
        %  This function finds the index (i.e. the path from the root to
        %  the patch recursively for all the leaves of this patch.
        %
        %Input:
        %   index    : the index from the root to this patch
        %Output:
        function findIndex(obj,index)
            
            for k=1:2
                index_c = [index k];
                if obj.children{k}.is_leaf
                    obj.children{k}.index = index_c;
                else
                    obj.children{k}.findIndex(index_c);
                end
            end
            
        end
        
        
        function vals = evalfGrid(obj,X,dim,order)
            [sum,dotprod] = obj.evalfGrid_recurse(X);
            vals = dotprod./sum;
        end
        
        function vals = evalf(obj,X,dim,order)
            [sum,dotprod] = obj.evalf_recurse(X);
            vals = dotprod./sum;
        end
        
        function [sum,dotprod] = evalf_recurse(obj,X)
            
            [num_pts,~] = size(X);
            
            dotprod = zeros(num_pts,1);
            
            sum = zeros(num_pts,1);
            
            %calculate values for the children
            for k=1:2
                ind = obj.children{k}.InDomain(X);
                
                
                if any(ind)
                    if ~obj.children{k}.is_leaf
                        [sumk,dotprodk] = obj.children{k}.evalf_recurse(X(ind,:));
                    else
                        sumk = obj.children{k}.evalfBump(X(ind,:));
                        dotprodk = sumk.*obj.children{k}.evalf(X(ind,:),1,0);
                    end
                        sum(ind) = sum(ind) + sumk;
                        dotprod(ind) = dotprod(ind) + dotprodk;
                end
            end
        end
        
        
        
        function [sum,dotprod] = evalfGrid_recurse(obj,X)

            grid_lengths = cellfun(@(x)length(x),X);
            
            sum = zeros(grid_lengths);
            
            dotprod = zeros(grid_lengths);
            
            %calculate values for the children
            for k=1:2
                
                [sub_grid,sub_ind] = obj.children{k}.IndDomainGrid(X);
                
                if all(cellfun(@any,sub_ind))
                    if ~obj.children{k}.is_leaf
                            [sumk,dotprodk] = obj.children{k}.evalfGrid_recurse(sub_grid);
                    else
                        sumk = obj.children{k}.evalfGridBump(sub_grid);
                        dotprodk = sumk.*obj.children{k}.evalfGrid(sub_grid);
                    end
                    
                    if obj.dim == 2
                            sum(sub_ind{1},sub_ind{2}) = sum(sub_ind{1},sub_ind{2}) + sumk;
                            dotprod(sub_ind{1},sub_ind{2}) = dotprod(sub_ind{1},sub_ind{2}) + dotprodk;
                    else
                            sum(sub_ind{1},sub_ind{2},sub_ind{3}) = sum(sub_ind{1},sub_ind{2},sub_ind{3}) + sumk;
                            dotprod(sub_ind{1},sub_ind{2},sub_ind{3}) = dotprod(sub_ind{1},sub_ind{2},sub_ind{3}) + dotprodk;
                    end
                end
            end
        end
        
        function M = interpSparseMatrixZone(obj,X)
            [ii,jj,zz] = obj.interpMatrixZone_vecs(X);
            M = sparse(ii,jj,zz,size(X,1),length(obj));
        end
            
        function [ii,jj,zz] = interpMatrixZoneGrid_vecs(obj,grid)
            
            ii = []; jj = []; zz = [];
            grid_lengths = cellfun(@(x)length(x),grid);
            
            [n,~] = size(grid_lengths);
            
            if n>1
                grid_lengths = grid_lengths';
            end
            
            sub_grids = cell(2,1);
            
            inds = cell(2,1);
            
            not_empty_child = true(2,1);
            
            sub_grids{1} = grid;
            sub_grids{2} = grid;
            
            midpoint = mean(obj.zone(obj.splitting_dim,:));
            
            inds{1} = grid{obj.splitting_dim} < midpoint;
            inds{2} = ~inds{1};
            
            for k=1:2
                sub_grids{k}{obj.splitting_dim} = sub_grids{k}{obj.splitting_dim}(inds{k});
                not_empty_child(k) = all(cellfun(@(x)~isempty(x),sub_grids{k}));
            end
            
            in_index = (1:prod(grid_lengths))';
            
            step = 0;
            
            for k=1:2
                if not_empty_child(k)
                    if obj.children{k}.is_leaf
                        [ii_k,jj_k,zz_k] = find(obj.children{k}.interpMatrixGrid(sub_grids{k}));
                        
                        if size(ii_k,2)~=1
                            ii_k = ii_k';
                        end
                        
                        if size(jj_k,2)~=1
                            jj_k = jj_k';
                        end
                        
                        if size(zz_k,2)~=1
                            zz_k = zz_k';
                        end
                        
                    else
                        [ii_k,jj_k,zz_k] = obj.children{k}.interpMatrixZoneGrid_vecs(sub_grids{k});
                    end
                    
                    if obj.dim==2
                        ind_k = false(grid_lengths(1),grid_lengths(2));
                        
                        if obj.splitting_dim==1
                            ind_k(inds{k},:) = true;
                        else
                            ind_k(:,inds{k}) = true;
                        end
                    elseif obj.dim==3
                        ind_k = false(grid_lengths(1),grid_lengths(2),grid_lengths(3));
                        
                        if obj.splitting_dim==1
                            ind_k(inds{k},:,:) = true;
                        elseif obj.splitting_dim==2
                            ind_k(:,inds{k},:) = true;
                        else
                            ind_k(:,:,inds{k}) = true;
                        end
                    end
                    
                    ind_k = ind_k(:);
                    
                    ind_k = in_index(ind_k);
                    
                    ii = [ii;ind_k(ii_k)];
                    jj = [jj;jj_k+step];
                    zz = [zz;zz_k];
                    step = step + obj.children{k}.length();
                end
                
            end
        end
        
        function vals = evalfZoneGrid(obj,X)
            grid_lengths = cellfun(@(x)length(x),X);
            
            [n,~] = size(grid_lengths);
            
            if n>1
                grid_lengths = grid_lengths';
            end
            
            %Put the order at the end; this makes things
            %easier to deal with in matlab
            vals = zeros(grid_lengths);
            
            sub_grids = cell(2,1);
            
            inds = cell(2,1);
            
            not_empty_child = true(2,1);
            
            sub_grids{1} = X;
            sub_grids{2} = X;
            
            midpoint = mean(obj.zone(obj.splitting_dim,:));
            
            inds{1} = X{obj.splitting_dim} < midpoint;
            inds{2} = ~inds{1};
            
            for k=1:2
                sub_grids{k}{obj.splitting_dim} = sub_grids{k}{obj.splitting_dim}(inds{k});
                not_empty_child(k) = all(cellfun(@(x)~isempty(x),sub_grids{k}));
            end
            
            
            %calculate values for the children
            for k=1:2
                if  not_empty_child(k)
                    if ~obj.children{k}.is_leaf
                        child_vals = obj.children{k}.evalfZoneGrid(sub_grids{k});
                    else
                        child_vals = obj.children{k}.evalfGrid(sub_grids{k});
                    end
                    
                    if obj.dim== 2
                        if obj.splitting_dim == 1
                            vals(inds{k},:) = child_vals;
                        else
                            vals(:,inds{k}) = child_vals;
                        end
                    elseif obj.dim == 3
                        if obj.splitting_dim == 1
                            vals(inds{k},:,:) = child_vals;
                        elseif obj.splitting_dim == 2
                            vals(:,inds{k},:) = child_vals;
                        else
                            vals(:,:,inds{k}) = child_vals;
                        end
                    end
                end
            end
        end
        
        function [ii,jj,zz] = interpMatrixZone_vecs(obj,X)
            
            [numpts,~] = size(X);
            
            ind = false(numpts,2);
            
            ii = [];
            jj = [];
            zz = [];
            
            mid_point = mean(obj.zone(obj.splitting_dim,:));
            
            ind(:,1) = X(:,obj.splitting_dim)<mid_point;
            ind(:,2) = ~ind(:,1);
            
            step = 0;
            
            %calculate values for the children
            for k=1:2
                if any(ind(:,k))
                    
                in_index = (1:numpts)';
                in_index = in_index(ind(:,k));
                
                    if obj.children{k}.is_leaf
                        M = obj.children{k}.interpMatrixPoints(X(ind(:,k),:));
                        [iik,jjk,zzk] = find(M);
                        
                        if(size(iik,2)~=1)
                            iik = iik.';
                            jjk = jjk.';
                            zzk = zzk.';
                        end
                        
                    else
                        [iik,jjk,zzk] = obj.children{k}.interpMatrixZone_vecs(X(ind(:,k),:));
                    end
                    
                    ii = [ii;in_index(iik)];
                    jj = [jj;jjk+step];
                    zz = [zz;zzk];
                end
                step = step+length(obj.children{k});
            end
        end
        
        function vals = evalfZone(obj,X)
            
            [numpts,~] = size(X);
            
            vals = zeros(numpts,1);
            
            ind = false(numpts,2);
            
            ind(:,1) = X(:,obj.splitting_dim)<mid_point;
            ind(:,2) = ~ind(:,1);
            
            child_vals = cell(2,1);
            
            %calculate values for the children
            for k=1:2
                if any(ind(:,k))
                    if obj.children{k}.is_leaf
                        child_vals{k} = obj.children{k}.evalf(X(ind(:,k),:));
                    else
                        child_vals{k} = obj.children{k}.evalfZone(X(ind(:,k),:));
                    end
                else
                    child_vals{k} = [];
                end
            end
            
            vals(ind(:,1)) = child_vals{1};
            vals(ind(:,2)) = child_vals{2};
        end
        
        
        function [LEAVES]= collectLeaves(obj,leaves)
            for k=1:2
                if obj.children{k}.is_leaf
                    leaves{length(leaves)+1} = obj.children{k};
                else
                    leaves = obj.children{k}.collectLeaves(leaves);
                end
            end
            LEAVES = leaves;
        end
        
        function IsGeometricallyRefined = IsGeometricallyRefined(obj)
            G1 = obj.children{1}.IsGeometricallyRefined();
            G2 = obj.children{2}.IsGeometricallyRefined();
            IsGeometricallyRefined = G1 & G2;
        end
        
        function plotdomain(obj)
            obj.children{1}.plotdomain();
            obj.children{2}.plotdomain();
        end
        
        function plotzone(obj)
            if(~isempty(obj.children{1}.zone))
                obj.children{1}.plotzone();
            end
            
            if(~isempty(obj.children{2}.zone))
                obj.children{2}.plotzone();
            end
        end
        

        
        function Coarsen(obj)
            for k=1:2
                obj.children{k}.Coarsen();
            end
        end
        
        function Refine(obj)
            
            
            for k=1:2
                obj.children{k}.Refine();
            end
        end
        
        function vals = Getvalues(obj)
            vals = [obj.children{1}.Getvalues();obj.children{2}.Getvalues()];
        end
        
        function split(obj)
            for k=1:2
                if obj.children{k}.is_leaf
                    lengths = diff(obj.children{k}.zone');
                    [~,split_dim] = max(lengths);
                    obj.children{k} = obj.children{k}.split(split_dim);
                else
                    obj.children{k}.split(Max,set_vals);
                end
            end
            obj.domain = [obj.children{1}.region(:,1) obj.children{2}.region(:,2)];
        end
        
        
        function str = toString(obj)
            str = strvcat(strcat('1',obj.children{1}.toString()),strcat('2',obj.children{2}.toString()));
        end
        
    end
    
end

