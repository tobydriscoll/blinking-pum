classdef (Abstract) Patch < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        domain
        boundaryIndex
        cheb_length
        is_leaf
        is_refined = false
        is_geometric_refined = false
        dim
        tol
        zone
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
            domain_ind = true(length(x),1);
            for i=1:obj.dim
                domain_ind = domain_ind & ...
                    ( x(:,i)>= obj.zone(i,1) & x(:,i)<= obj.zone(i,2));
            end
        end
        
                %Right now just assume the domain is just the square.
        function [sub_grid,split_ind] = InDomainGrid(obj,x,split_dim)
            sub_grid = cell(1,obj.dim);
            for i=1:obj.dim
                ind = x{i}>=obj.domain(i,1) & x{i}<=obj.domain(i,2);
                
                if i==split_dim
                    split_ind = ind;
                end
                
                sub_grid{i} = x{i}(ind);
            end
        end 
    end
end

