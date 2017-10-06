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
        
        function is_refined = PUsplit(obj)
            
            is_refined = true;
            is_geometric_refined = true;
            
            for k=1:2
                
                if ~obj.children{k}.is_refined
                    if obj.children{k}.is_leaf
                        obj.children{k} = obj.children{k}.splitleaf();
                        is_refined = is_refined && obj.children{k}.is_refined;
                        is_geometric_refined = is_geometric_refined && obj.children{k}.is_geometric_refined;
                    else
                        temp = obj.children{k}.PUsplit();
                        is_refined = is_refined && temp;
                        is_geometric_refined = is_geometric_refined && obj.children{k}.is_geometric_refined;
                    end
                end
            end
            
            obj.region = [obj.children{1}.region(:,1) obj.children{2}.region(:,2)];
            obj.domain = obj.region;
            
            obj.cheb_length = obj.children{1}.cheb_length+obj.children{2}.cheb_length;
            obj.overlap_in = [obj.children{2}.region(obj.splitting_dim,1), obj.children{1}.region(obj.splitting_dim,2)];
            
            obj.is_refined = is_refined;
            obj.is_geometric_refined = is_geometric_refined;
        end
        
        function sample(obj,f)
            if ~isnumeric(f)
                obj.children{1}.sample(f);
                obj.children{2}.sample(f);
            else
                obj.children{1}.sample(f(1:length(obj.children{1})));
                obj.children{2}.sample(f(length(obj.children{1})+1:end));
            end
        end
        
        
        %  [WEIGHTSVALS]= evalweights(obj,index,X,dim,order)
        %Input:
        %   index    : the index for the leaf of the weight to be evaluated
        %   x        : the set of points in each dimension.
        %   dim      : the dimension the derivative will be taken in
        %   order    : order of the derivative
        %Output:
        %   WEIGHTSVALS :  cell array of the weightvalues in each
        %                  dimension, including the dimension
        function [WEIGHTSVALS]= evalweights(obj,index,X,dim,order)
            %Our weights are seperable, so we compute the weights
            %for the functions in each dimension seperately.
            if obj.children{index(1)}.is_leaf
                
                %as a base case, we set the weight to one.
                for i=1:obj.dim
                    WEIGHTSVALS{i} = zeros(length(X{i}),order+1);
                    WEIGHTSVALS{i}(:,1) = ones(length(X{i}),1);
                end
                
                WEIGHTSVALS{obj.splitting_dim}(:,1) = obj.weights.evalf(X{obj.splitting_dim},obj.overlap_in,index(1),1,0);
            else
                CHILDWEIGHT = obj.children{index(1)}.evalweights(index(2:end),X,dim,order);
                WEIGHTSVALS = CHILDWEIGHT;
                
                WEIGHTSVALS{obj.splitting_dim} = zeros(length(X{obj.splitting_dim}),order+1);
                
                %Use the product role to compute the derivative of the
                %weight.
                for k=0:order
                    for j=0:k
                        if (k-j)==0 || obj.splitting_dim == dim
                            WEIGHTSVALS{obj.splitting_dim}(:,k+1) =  WEIGHTSVALS{obj.splitting_dim}(:,k+1) + ...
                                nchoosek(k,j).*obj.weights.evalf(X{obj.splitting_dim},obj.overlap_in,index(1),1,k-j).*CHILDWEIGHT{obj.splitting_dim}(:,j+1);
                        end
                    end
                end
            end
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
                        [sumk,dotprodk] = obj.children{k}.evalf(X(ind,:),dim,order);
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
                        dotprodk = sumk.*obj.children{k}.evalfGrid(sub_grid,1,0);
                    end
                    
                    if obj.dim ==2
                            sum(sub_ind{1},sub_ind{2}) = sum(sub_ind{1},sub_ind{2}) + sumk;
                            dotprod(sub_ind{1},sub_ind{2}) = dotprod(sub_ind{1},sub_ind{2}) + dotprodk;
                    else
                            sum(sub_ind{1},sub_ind{2},sub_ind{3}) = sum(sub_ind{1},sub_ind{2},sub_ind{3}) + sumk;
                            dotprod(sub_ind{1},sub_ind{2},sub_ind{3}) = dotprod(sub_ind{1},sub_ind{2},sub_ind{3}) + dotprodk;
                    end
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
            lengths = size(vals);
            
            sub_grids = cell(2,1);
            child_vals = cell(2,1);
            
            inds = cell(2,1);
            
            not_empty_child = true(2,1);
            
            sub_grids{1} = X;
            sub_grids{2} = X;
            
            midpoint = mean(obj.zone(obj.splitting_dim,:));
            
            inds{1} = X{obj.splitting_dim} < midpoint;
            inds{2} = X{obj.splitting_dim} >= midpoint;
            
            for k=1:2
                sub_grids{k}{obj.splitting_dim} = sub_grids{k}{obj.splitting_dim}(inds{k});
                not_empty_child(k) = all(cellfun(@(x)~isempty(x),sub_grids{k}));
            end
            
            
            %calculate values for the children
            for k=1:2
                if  not_empty_child(k)
                    if ~obj.children{k}.is_leaf
                        child_vals{k} = obj.children{k}.evalfZoneGrid(sub_grids{k});
                    else
                        child_vals{k} = obj.children{k}.evalfGrid(sub_grids{k},1,0);
                    end
                end
            end
            
            %Go ahead and shift the splitting diminsion
            %Here I ciculate the dimensions before the order.
            dim_permute = circshift(1:obj.dim,[0 -obj.splitting_dim+1]);
            
            vals = permute(vals,dim_permute);
            lengths = size(vals);
            vals = reshape(vals,lengths(1),prod(lengths(2:end)));
            
            for k=1:2
                child_vals{k} = permute(child_vals{k},[dim_permute obj.dim+1]);
                c_lengths = size(child_vals{k});
                child_vals{k} = reshape(child_vals{k},c_lengths(1),prod(c_lengths(2:end)));
            end            
            
            for k=1:2
                if not_empty_child(k)
                    vals(inds{k},:) = child_vals{k};
                end
            end
            
            vals = reshape(vals,lengths);
            vals = ipermute(vals,dim_permute);
        end
        
        function vals = evalfZone(obj,X)
            
            [numpts,~] = size(X);
            
            vals = zeros(numpts,1);
            
            ind = false(numpts,2);
            
            %Figure out indicies of the points in the left and right
            %domains
            for k=1:2
                ind(:,k) = obj.children{k}.InZone(X);
            end
            
            child_vals = cell(2,1);
            
            %calculate values for the children
            for k=1:2
                if any(ind(:,k))
                    if obj.children{k}.is_leaf
                        child_vals{k} = obj.children{k}.evalf(X(ind(:,k),:),1,0);
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
        
        function ResolveChebWeights(obj,weights)
            
            w = weights{obj.splitting_dim};
            
            for k=1:2
                
                domain = obj.domain(obj.splitting_dim,:);
                subdomain = obj.children{k}.domain(obj.splitting_dim,:);
                h = @(x) 2/diff(domain)*x-sum(domain)/diff(domain);
                
                wk = obj.weights.chebweights(:,k);
                wk = wk(chebfun(h,domain));
                
                weights{obj.splitting_dim} = chebfun(@(x)w(x).*wk(x),subdomain);
                
                if obj.children{k}.is_leaf
                    obj.children{k}.chebweights = weights;
                else
                    obj.children{k}.ResolveChebWeights(weights);
                end
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
        
        function split(obj,overlap)
            for k=1:2
                if obj.children{k}.is_leaf
                    lengths = diff(obj.children{k}.domain');
                    [~,split_dim] = max(lengths);
                    obj.children{k} = obj.children{k}.split(split_dim);
                else
                    obj.children{k}.split(overlap);
                end
            end
        end
        
        
        function str = toString(obj)
            str = strvcat(strcat('1',obj.children{1}.toString()),strcat('2',obj.children{2}.toString()));
        end
    end
    
end

