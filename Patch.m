classdef (Abstract) Patch < handle & matlab.mixin.Copyable
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        domain
        zone
        boundaryIndex
        cheb_length
        Root
        is_leaf
        is_refined = false
        is_geometric_refined = false
        dim
        tol
        overlap = 0.08;
        outerbox
        split_flag %array indicating if we will split along a dimension
    end
    
    methods (Abstract)
        points = points(obj)
        ef = evalf(obj,x,dim,order)
        ef = evalfGrid(obj,x,dim,order)
        ln = length(obj)
        sample(obj,f)
        plotdomain(obj)
    end
    
    methods
        %Right now just assume the domain is just the square.
        function domain_ind = InDomain(obj,x)
            [num_pts,~] = size(x);
            domain_ind = true(num_pts,1);
            for i=1:obj.dim
                domain_ind = domain_ind & ...
                    ( x(:,i)>= obj.domain(i,1) & x(:,i)<= obj.domain(i,2));
            end
        end
        
        %Right now just assume the domain is just the square.
        function domain_ind = InZone(obj,x)
            domain_ind = true(size(x,1),1);
            for i=1:obj.dim
                domain_ind = domain_ind & ...
                    ( x(:,i)>= obj.zone(i,1) & x(:,i)<= obj.zone(i,2));
            end
        end
        
        %Right now just assume the domain is just the square.
        function [sub_grid,split_ind] = InDomainGrid(obj,x,split_dim)
            sub_grid = cell(1,obj.dim);
            split_grid = cell(1,obj.dim);
            
            for i=1:obj.dim
                ind = x{i}>=obj.domain(i,1) & x{i}<=obj.domain(i,2);
                
                if i==split_dim
                    split_ind = ind;
                end
                
                sub_grid{i} = x{i}(ind);
            end
        end
        
        function [sub_grid,split_ind] = IndDomainGrid(obj,x)
            sub_grid = cell(1,obj.dim);
            split_ind = cell(1,obj.dim);
            
            for i=1:obj.dim
                ind = x{i}>=obj.domain(i,1) & x{i}<=obj.domain(i,2);
                
                split_ind{i} = ind;
                
                sub_grid{i} = x{i}(ind);
            end
        end
        
        function show(obj,level)
            % Make a pretty graph showing the domains (2D)
            assert(obj.dim==2,'Must be 2-D')
            if nargin==1  % user call
                level = 0;
                newplot
            end
            if obj.is_leaf
                x = obj.domain(1,[1 1 2 2]);
                y = obj.domain(2,[1 2 2 1]);
                z = level*ones(1,4);
                patch(x,y,z,rand(1,3),'facealpha',0.2,'edgecolor','none')
            else
                show(obj.children{1},level+1);
                show(obj.children{2},level+1);
            end
        end
        
        function clean(obj)
            if obj.is_leaf
                obj.split_flag = false(size(obj.split_flag));
            else
                
                clean(obj.children{1});
                clean(obj.children{2});
                obj.split_flag = obj.children{1}.split_flag | obj.children{2}.split_flag;
                obj.split_flag(obj.splitting_dim) = true;
                
                for i=1:obj.dim
                    obj.domain(i,:) = [min(obj.children{1}.domain(i,1),min(obj.children{2}.domain(i,1))) max(obj.children{1}.domain(i,2),min(obj.children{2}.domain(i,2)))];
                end
                
                obj.cheb_length = obj.children{1}.cheb_length()+obj.children{2}.cheb_length();
            end
        end
        
        function fun_obj = refine(obj,f,grid_opt)
            
            if nargin<3
                grid_opt = false;
            end
            
            fun_obj = obj;
            
            while ~fun_obj.is_refined
                
                Max = fun_obj.sample(f,grid_opt);
                
                if fun_obj.is_leaf
                    fun_obj = obj.splitleaf(Max);
                else
                    fun_obj.PUsplit(Max);
                end
                
            end
            
        end
        
        function isGeoRefined = IsGeometricallyRefined(obj)
            
            if obj.is_leaf
                isGeoRefined = IsLeafGeometricallyRefined(obj);
                obj.is_geometric_refined = isGeoRefined;
            else
                isGeoRefined = true;
                for k=1:2
                    isGeoRefined = isGeoRefined && IsGeometricallyRefined(obj.children{k});
                end
            end
            
        end
        
        
        
        
    end
end

