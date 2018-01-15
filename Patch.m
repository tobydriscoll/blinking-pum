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
                    fun_obj = obj.splitleaf(Max);
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
                
                T_add = PUPatch(domain_add,obj.zone,cheb_length_add,children,obj.splitting_dim,obj.index);
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
        
        function T_sub = subtract(obj,T_2,sub_f)
            
            if obj.is_leaf
                
                T_sub = copy(T_2);
                
                if ~T_sub.is_leaf
                    T_sub.findIndex([]);
                    leafArray = T_sub.collectLeaves({});
                else
                    leafArray = {T_sub};
                end
                
                for i=1:length(leafArray)
                    if isequal(obj.domain,leafArray{i}.domain) && isequal(obj.degs,leafArray{i}.degs)
                        leafArray{i}.values = obj.values-leafArray{i}.values;
                    else
                        leafArray{i}.values = obj.evalfGrid(leafArray{i}.leafGrids())-leafArray{i}.values;
                    end
                end
                
            elseif T_2.is_leaf
                T_sub = copy(obj);
                
                if ~T_sub.is_leaf
                    T_sub.findIndex([]);
                    leafArray = T_sub.collectLeaves({});
                else
                    leafArray = {T_sub};
                end
                
                for i=1:length(leafArray)
                    if isequal(T_2.domain,leafArray{i}.domain) && isequal(T_2.degs,leafArray{i}.degs)
                        leafArray{i}.values = leafArray{i}.values-T_2.values;
                    else
                        leafArray{i}.values = leafArray{i}.values-T_2.evalfGrid(leafArray{i}.leafGrids());
                    end
                end
                
            elseif obj.splitting_dim == T_2.splitting_dim
                
                children{1} = subtract(obj.children{1},T_2.children{1},sub_f);
                children{2} = subtract(obj.children{2},T_2.children{2},sub_f);
                
                cheb_length_add = length(obj.children{1})+length(obj.children{2});
                domain_add = [obj.children{1}.domain(:,1) obj.children{2}.domain(:,2)];
                
                T_sub = PUPatch(domain_add,obj.zone,cheb_length_add,children,obj.splitting_dim,obj.index);
            else
                leafArray = obj.collectLeaves({});
                num_patch1 = length(leafArray);
                
                leafArray = obj.collectLeaves({});
                num_patch2 = length(leafArray);
                
                if num_patch1>num_patch2
                    T_sub = copy(obj);
                else
                    T_sub = copy(T_2);
                end
                
                T_sub.reset();
                T_sub = T_sub.refine(sub_f,true);
            end
            
        end
        
        function T_mult = multiply(obj,T_2,mult_f)
            
            if obj.is_leaf
                
                T_mult = copy(T_2);
                
                if ~T_mult.is_leaf
                    T_mult.findIndex([]);
                    leafArray = T_mult.collectLeaves({});
                    T_2Array = T_2.collectLeaves({});
                else
                    leafArray = {T_mult};
                    T_2Array = {T_2};
                end
                
                T_mult.reset();
                
                for i=1:length(leafArray)
                    leafArray{i} = refine(leafArray{i},@(x)obj.evalfGrid(x).*T_2Array{i}.evalfGrid(x),true);
                end
                
            elseif T_2.is_leaf
                T_mult = copy(obj);
                
                if ~T_mult.is_leaf
                    T_Array = obj.collectLeaves({});
                    T_mult.findIndex([]);
                    leafArray = T_mult.collectLeaves({});
                else
                    T_Array = {obj};
                    leafArray = {T_mult};
                end
                
                T_mult.reset();
                
                for i=1:length(leafArray)
                    leafArray{i} = refine(leafArray{i},@(x)T_2.evalfGrid(x).*T_Array{i}.evalfGrid(x),true);
                end
                
            elseif obj.splitting_dim == T_2.splitting_dim
                
                children{1} = multiply(obj.children{1},T_2.children{1},mult_f);
                children{2} = multiply(obj.children{2},T_2.children{2},mult_f);
                
                cheb_length_add = length(obj.children{1})+length(obj.children{2});
                domain_add = [obj.children{1}.domain(:,1) obj.children{2}.domain(:,2)];
                
                T_mult = PUPatch(domain_add,obj.zone,cheb_length_add,children,obj.splitting_dim,obj.index);
            else
                leafArray = obj.collectLeaves({});
                num_patch1 = length(leafArray);
                
                leafArray = obj.collectLeaves({});
                num_patch2 = length(leafArray);
                
                if num_patch1>num_patch2
                    T_mult = copy(obj);
                else
                    T_mult = copy(T_2);
                end
                
                T_mult.reset();
                T_mult = T_mult.refine(mult_f,true);
            end
        end
        
        function T_divide = divide(obj,T_2,divide_f)
            
            if obj.is_leaf
                
                T_divide = copy(T_2);
                
                if ~T_divide.is_leaf
                    T_divide.findIndex([]);
                    leafArray = T_divide.collectLeaves({});
                    T_2Array = T_2.collectLeaves({});
                else
                    leafArray = {T_divide};
                    T_2Array = {T_2};
                end
                
                T_divide.reset();
                
                for i=1:length(leafArray)
                    leafArray{i} = refine(leafArray{i},@(x)obj.evalfGrid(x)./T_2Array{i}.evalfGrid(x),true);
                end
                
            elseif T_2.is_leaf
                T_divide = copy(obj);
                
                if ~T_divide.is_leaf
                    T_Array = obj.collectLeaves({});
                    T_divide.findIndex([]);
                    leafArray = T_divide.collectLeaves({});
                else
                    T_Array = {obj};
                    leafArray = {T_divide};
                end
                
                T_divide.reset();
                
                for i=1:length(leafArray)
                    leafArray{i} = refine(leafArray{i},@(x)T_Array{i}.evalfGrid(x)./T_2.evalfGrid(x),true);
                end
                
            elseif obj.splitting_dim == T_2.splitting_dim
                
                children{1} = divide(obj.children{1},T_2.children{1},divide_f);
                children{2} = divide(obj.children{2},T_2.children{2},divide_f);
                
                cheb_length_add = length(obj.children{1})+length(obj.children{2});
                domain_add = [obj.children{1}.domain(:,1) obj.children{2}.domain(:,2)];
                
                T_divide = PUPatch(domain_add,obj.zone,cheb_length_add,children,obj.splitting_dim,obj.index);
            else
                leafArray = obj.collectLeaves({});
                num_patch1 = length(leafArray);
                
                leafArray = obj.collectLeaves({});
                num_patch2 = length(leafArray);
                
                if num_patch1>num_patch2
                    T_divide = copy(obj);
                else
                    T_divide = copy(T_2);
                end
                
                T_divide.reset();
                T_divide = T_divide.refine(divide_f,true);
            end
        end
        
        function T_power = power(obj,p)
            T_power = copy(obj);
            
                T_power.reset();
                
                if ~T_power.is_leaf
                    T_Array = obj.collectLeaves({});
                    T_power.findIndex([]);
                    leafArray = T_power.collectLeaves({});
                else
                    T_Array = {obj};
                    leafArray = {T_power};
                end
            
            for i=1:length(T_Array)
                leafArray{i} = refine(leafArray{i},@(x)T_Array{i}.evalfGrid(x).^p,true);
            end
        end
        
        
    end
end

