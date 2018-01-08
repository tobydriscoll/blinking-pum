classdef (Abstract) Patch < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        domain
        zone
        boundaryIndex
        cheb_length
        is_leaf
        is_refined = false
        is_geometric_refined = false
        dim
        tol
        overlap = 0.08;
        outerbox
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
        
    end
end

