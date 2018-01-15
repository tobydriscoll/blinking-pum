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
        overlap = 0.1;
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
        
        function fun_obj = refine(obj,f,grid_opt)
            
            if nargin<3
                grid_opt = false;
            end
            
            fun_obj = obj;
            
            while ~fun_obj.is_refined
                
                Max = fun_obj.sample(f,grid_opt);
                
                if fun_obj.is_leaf
                    fun_obj = obj.ChebRoot.splitleaf(Max);
                else
                    fun_obj.PUsplit(Max);
                end
                
            end
            
        end
        
        function T_add = add(obj,T_2,add_f)
            
            if obj.is_leaf
                
                T_add = copy(T_2);
                
                if ~T_add.is_leaf
                    T_add.findIndex([]);
                    leafArray = T_add.collectLeaves({});
                else
                    leafArray = {T_add};
                end
                
                for i=1:length(leafArray)
                    if isequal(obj.domain,leafArray{i}.domain) && isequal(obj.degs,leafArray{i}.degs)
                        leafArray{i}.values = leafArray{i}.values+obj.values;
                    else
                        leafArray{i}.values = leafArray{i}.values+obj.evalfGrid(leafArray{i}.leafGrids());
                    end
                end
                
            elseif T_2.is_leaf
                T_add = copy(obj);
                
                if ~T_add.is_leaf
                    T_add.findIndex([]);
                    leafArray = T_add.collectLeaves({});
                else
                    leafArray = {T_add};
                end
                
                for i=1:length(leafArray)
                    if isequal(T_2.domain,leafArray{i}.domain) && isequal(T_2.degs,leafArray{i}.degs)
                        leafArray{i}.values = leafArray{i}.values+T_2.values;
                    else
                        leafArray{i}.values = leafArray{i}.values+T_2.evalfGrid(leafArray{i}.leafGrids());
                    end
                end
                
            elseif obj.splitting_dim == T_2.splitting_dim
                
                children{1} = add(obj.children{1},T_2.children{1},add_f);
                children{2} = add(obj.children{2},T_2.children{2},add_f);
                
                cheb_length_add = length(obj.children{1})+length(obj.children{2});
                domain_add = [obj.children{1}.domain(:,1) obj.children{2}.domain(:,2)];
                
                T_add = PUPatch(domain_add,obj.zone,[],cheb_length_add,children,obj.splitting_dim,obj.index);
            else
                leafArray = obj.collectLeaves({});
                num_patch1 = length(leafArray);
                
                leafArray = obj.collectLeaves({});
                num_patch2 = length(leafArray);
                
                if num_patch1>num_patch2
                    T_add = copy(obj);
                else
                    T_add = copy(T_2);
                end
                
                T_add.reset();
                T_add = T_add.refine(add_f,true);
            end
            
        end
        
        
    end
end

