classdef PUPatch<Patch
    % This class is represents a tree with two Patch Children.
    % Children can be any Patch class, including PUPatch and 
    % LeafPatch.
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
        
        % Construct for the PUPatch
        %
        %        Input:
        %        domain: (dim x 2) array indiciating array for the domain.
        %          zone: (dim x 2) array indiciating array for the zone.
        %    overlap_in: overlap interval along splitting dimension.
        %   cheb_length: total number of interpolating points.
        %      children: (2 x 1) cell array of Patch objects;
        % splitting_dim: dimension the patch is split in.
        %         index: array indicating patch from root to patch (1 left,
        %         right left).
        function obj = PUPatch(domain,zone,overlap_in,cheb_length,children,splitting_dim,index)
            obj.outerbox = children{1}.outerbox;
            obj.zone = zone;
            obj.domain = domain;
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
        
        % Returns the length of the patch
        %   Output:
        %       ln: total number of interpolating points
        function ln = length(obj)
            ln = obj.cheb_length;
        end
        
        % Returns cell array of grids on leaves
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
        
        % Returns vector of interpolating points for the leaf.
        function pts = points(obj)
            pts = [obj.children{1}.points();obj.children{2}.points()];
        end
        
        % Will split the children of the patch if they are unrefined.
        %
        %     Output: 
        % is_refined: returns if the patch is refined.
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
            
            obj.domain = [obj.children{1}.domain(:,1) obj.children{2}.domain(:,2)];
            
            obj.cheb_length = obj.children{1}.cheb_length+obj.children{2}.cheb_length;
            obj.overlap_in = [obj.children{2}.domain(obj.splitting_dim,1), obj.children{1}.domain(obj.splitting_dim,2)];
            
            obj.is_refined = is_refined;
            obj.is_geometric_refined = is_geometric_refined;
        end
        
        % This method samples the leaves of the patch.
        %
        %     Input:
        %         f: vector of values for the interpolating values on the 
        %            leaves depth first, or an anonymous function.
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
        
        %  evalfGrid(obj,X)
        %  Evaluates the approximation on a grid.
        %
        %Input:
        %    X: cell array of grid values.
        function vals = evalfGrid(obj,X)
            [sum,dotprod] = obj.evalfGrid_recurse(X);
            vals = dotprod./sum;
        end
        
        %  evalf(obj,X)
        %  Evaluates the approximation on a list of points.
        %
        %Input:
        %    X: list of points.
        function vals = evalf(obj,X)
            [sum,dotprod] = obj.evalf_recurse(X);
            vals = dotprod./sum;
        end
        
        %  evalf_recurse(obj,X)
        %  Evaluates the values need recursively for the approximation on a 
        %  list of points X.
        %
        %Input:
        %    X: list of points.
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
                        dotprodk = sumk.*obj.children{k}.evalf(X(ind,:));
                    end
                        sum(ind) = sum(ind) + sumk;
                        dotprod(ind) = dotprod(ind) + dotprodk;
                end
            end
        end
        
        %  evalfGrid_recurse(obj,X)
        %  Evaluates the values need recursively for the approximation on a 
        %  grid X.
        %
        %Input:
        %    X: cell array of grid values.
        function [sum,dotprod] = evalfGrid_recurse(obj,X)
            
            grid_lengths = cellfun(@(x)size(x,1),X);
            
            sum = zeros(grid_lengths);
            
            dotprod = zeros(grid_lengths);
            
            for k=1:2
                [sub_grid{k},sub_ind{k}] = obj.children{k}.IndDomainGrid(X);
            end
            
            for i=1:obj.dim
                com_ind{i} = sub_ind{1}{i} & sub_ind{2}{i};
                
                for k=1:2
                    com_sub_ind{k}{i} = com_ind{i}(sub_ind{k}{i});
                end
            end
            
            
            %             any1 = all(cellfun(@any,sub_ind{1}));
            %             any2 = all(cellfun(@any,sub_ind{2}));
            %
            %             if any1
            %                 if ~obj.children{1}.is_leaf
            %                     [sum1,dotprod1] = obj.children{1}.evalfGrid_recurse(sub_grid{1});
            %                 else
            %                     sum1 = obj.children{1}.evalfGridBump(sub_grid{1});
            %                     dotprod1 = sum1.*obj.children{1}.evalfGrid(sub_grid{1});
            %                 end
            %             end
            %
            %             if any2
            %
            %                 if ~obj.children{2}.is_leaf
            %                     [sum2,dotprod2] = obj.children{2}.evalfGrid_recurse(sub_grid{2});
            %                 else
            %                     sum2 = obj.children{2}.evalfGridBump(sub_grid{2});
            %                     dotprod2 = sum2.*obj.children{2}.evalfGrid(sub_grid{2});
            %                 end
            %             end
            %
            %             if any1
            %
            %                 subses1 = repmat({':'}, [obj.dim]);
            %                 subses2 = repmat({':'}, [obj.dim]);
            %
            %                 for i=1:obj.dim
            %                     subses1{i} = sub_ind{1}{i} & ~ com_ind{i};
            %                     subses2{i} = ~com_sub_ind{1}{i};
            %                 end
            %
            %                 sum(subses1{:}) = sum1(subses2{:});
            %                 dotprod(subses1{:}) = dotprod1(subses2{:});
            %             end
            %
            %             if any1 || any2
            %
            %                 subses1 = repmat({':'}, [obj.dim]);
            %                 subses2 = repmat({':'}, [obj.dim]);
            %                 subses3 = repmat({':'}, [obj.dim]);
            %
            %                 for i=1:obj.dim
            %                     subses1{i} = com_ind{i};
            %                     subses2{i} = com_sub_ind{1}{i};
            %                     subses3{i} = com_sub_ind{2}{i};
            %                 end
            %
            %                 sum(subses1{:}) = sum1(subses2{:})+sum2(subses3{:});
            %                 dotprod(subses1{:}) = dotprod1(subses2{:})+dotprod2(subses3{:});
            %             end
            %
            %             if any2
            %
            %                 subses1 = repmat({':'}, [obj.dim]);
            %                 subses2 = repmat({':'}, [obj.dim]);
            %
            %                 for i=1:obj.dim
            %                     subses1{i} = sub_ind{2}{i} & ~ com_ind{i};
            %                     subses2{i} = ~com_sub_ind{2}{i};
            %                 end
            %
            %                 sum(subses1{:}) = sum2(subses2{:});
            %                 dotprod(subses1{:}) = dotprod2(subses2{:});
            %
            %             end
            
            %            calculate values for the children
            
            for k=1:2
                if all(cellfun(@any,sub_ind{k}))
                    if ~obj.children{k}.is_leaf
                        [sumk,dotprodk] = obj.children{k}.evalfGrid_recurse(sub_grid{k});
                    else
                        sumk = obj.children{k}.evalfGridBump(sub_grid{k});
                        dotprodk = sumk.*obj.children{k}.evalfGrid(sub_grid{k});
                    end
                    
                    if obj.dim == 2
                        
                        if k==1
                            sum(sub_ind{k}{1},sub_ind{k}{2}) = sumk;
                            dotprod(sub_ind{k}{1},sub_ind{k}{2}) = dotprodk;
                            
                        else
                            
                            sum(sub_ind{k}{1} & ~ com_ind{1},sub_ind{k}{2} & ~ com_ind{2}) = sumk(~com_sub_ind{k}{1},~com_sub_ind{k}{2});
                            dotprod(sub_ind{k}{1} & ~ com_ind{1},sub_ind{k}{2} & ~ com_ind{2}) = dotprodk(~com_sub_ind{k}{1},~com_sub_ind{k}{2});
                            
                            sum(com_ind{1},com_ind{2}) = sum(com_ind{1},com_ind{2})+sumk(com_sub_ind{k}{1},com_sub_ind{k}{2});
                            dotprod(com_ind{1},com_ind{2}) = dotprod(com_ind{1},com_ind{2})+dotprodk(com_sub_ind{k}{1},com_sub_ind{k}{2});
                        end
                        
                    else
                        
                        if k==1
                            sum(sub_ind{k}{1},sub_ind{k}{2},sub_ind{k}{3}) = sumk;
                            dotprod(sub_ind{k}{1},sub_ind{k}{2},sub_ind{k}{3}) = dotprodk;
                            
                        else
                            
                            sum(sub_ind{k}{1} & ~ com_ind{1},sub_ind{k}{2} & ~ com_ind{2},sub_ind{k}{3} & ~ com_ind{3}) = sumk(~com_sub_ind{k}{1},~com_sub_ind{k}{2},~com_sub_ind{k}{3});
                            dotprod(sub_ind{k}{1} & ~ com_ind{1},sub_ind{k}{2} & ~ com_ind{2},sub_ind{k}{3} & ~ com_ind{3}) = sumk(~com_sub_ind{k}{1},~com_sub_ind{k}{2},~com_sub_ind{k}{3});
                            
                            sum(com_ind{1},com_ind{2},com_ind{3}) = sum(com_ind{1},com_ind{2},com_ind{3})+sumk(com_sub_ind{k}{1},com_sub_ind{k}{2},com_sub_ind{k}{3});
                            dotprod(com_ind{1},com_ind{2},com_ind{3}) = dotprod(com_ind{1},com_ind{2},com_ind{3})+dotprodk(com_sub_ind{k}{1},com_sub_ind{k}{2},com_sub_ind{k}{3});
                        end
                        
                    end
                end
            end
        end
        
        %  interpSparseMatrixZone(obj,X)
        %  Creates a sparse interpolating matrix for a list of points X.
        %
        % Output:
        %      M: sparse interpolating matrix.
        function M = interpSparseMatrixZone(obj,X)
            [ii,jj,zz] = obj.interpMatrixZone_vecs(X);
            M = sparse(ii,jj,zz,size(X,1),length(obj));
        end
        
        %  interpSparseMatrixZoneGrid(obj,X)
        %  Creates a sparse interpolating matrix along the zones
        %  for a grid X.
        %
        % Output:
        %      M: sparse interpolating matrix.
        function M = interpSparseMatrixZoneGrid(obj,X)
            [ii,jj,zz] = obj.interpMatrixZoneGrid_vecs(X);
            grid_lengths = cellfun(@(x)length(x),X);
            M = sparse(ii,jj,zz,prod(grid_lengths),length(obj));
        end
        
        %  interpSparseMatrixZone_vecs(obj,X)
        %  Creates a sparse interpolating matrix along the zones
        %  for a list of points X.
        %
        %     Input:
        %      grid: cell array of grid values.
        %
        %    Output:
        %  ii,jj,kk: vectors used to construct the sparse matrix.
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
                    
                end
                step = step + obj.children{k}.length();
            end
        end
        
        %  evalfZoneGrid(obj,X)
        %  Evaluates a grid along the zones for a grid X.
        %
        %     Input:
        %      grid: cell array of grid values.
        %
        %    Output:
        %      vals: values of approximation on grid along the zones.
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
        
        %  interpMatrixZone_vecs(obj,X)
        %  Creates a sparse interpolating matrix for a list of points X.
        %
        %     Input:
        %      grid: cell array of grid values.
        %
        %    Output:
        %  ii,jj,kk: vectors used to construct the sparse matrix.
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
        
        %  evalfZone(obj,X)
        %  Evaluates the approximation along the zones for a list of
        %  points.
        %
        %     Input:
        %         X: list of points.
        %
        %    Output:
        %      vals: values of approximation along the zones 
        %            for the list of points.
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
        
        %  collectLeaves(obj,leaves)
        %  Recursive function that collects leaves into cell array.
        %
        %     Input:
        %    leaves: current list of leaves
        %
        %    Output:
        %    LEAVES: leaves list with children of patch added.
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
        
        % Plots the domains of the children.
        function plotdomain(obj)
            obj.children{1}.plotdomain();
            obj.children{2}.plotdomain();
        end
        
        % Plots the zones of the children.
        function plotzone(obj)
            if(~isempty(obj.children{1}.zone))
                obj.children{1}.plotzone();
            end
            
            if(~isempty(obj.children{2}.zone))
                obj.children{2}.plotzone();
            end
        end
        
        %Coarsens the leaves of the patch.
        function Coarsen(obj)
            for k=1:2
                obj.children{k}.Coarsen();
            end
            obj.cheb_length = length(obj.children{1})+length(obj.children{2});
        end
        
        %Refines the leaves of the patch.
        function Refine(obj)
            for k=1:2
                obj.children{k}.Refine();
            end
            obj.cheb_length = length(obj.children{1})+length(obj.children{2});
        end
        
        %Method returns vector of values of interpolating points of patch.
        function vals = Getvalues(obj)
            vals = [obj.children{1}.Getvalues();obj.children{2}.Getvalues()];
        end
        
        %Recursive method that splits the children of a patch along a given
        %dimension.
        %
        %       Input:
        %   split_dim: splitting dimension
        %    set_vals: indicator if new children will have values 
        %              interpolated from the parent.
        function split(obj,split_dim,set_vals)
            for k=1:2
                if obj.children{k}.is_leaf
                    obj.children{k} = obj.children{k}.split(split_dim,set_vals);
                else
                    obj.children{k}.split(split_dim,set_vals);
                end
            end
            obj.cheb_length = length(obj.children{1})+length(obj.children{2});
            obj.domain = [obj.children{1}.domain(:,1) obj.children{2}.domain(:,2)];
        end
        
        % String function for method.
        function str = toString(obj)
            str = strvcat(strcat('1',obj.children{1}.toString()),strcat('2',obj.children{2}.toString()));
        end
        
    end
    
end

