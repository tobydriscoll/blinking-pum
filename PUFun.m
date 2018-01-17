classdef PUFun < handle & matlab.mixin.Copyable
    
    properties
        ChebRoot
        TreeGrid
        leafArray
        Errs
        Nums
        deg_in
        tol
        domain
    end
    
    methods
        
        function obj = PUFun(domain,deg_in,f,tol,grid_opt,ChebRoot)
            
            obj.domain = domain;
            obj.deg_in = deg_in;
            [dim,~] = size(domain);
            
            if nargin < 5
                obj.tol = 1e-12;
                grid_opt = false;
                obj.ChebRoot = ChebPatch(domain,domain,domain,deg_in,true(1,dim),obj.tol);
            elseif nargin < 6
                obj.tol = tol;
                obj.ChebRoot = ChebPatch(domain,domain,domain,deg_in,true(1,dim),obj.tol);
            else
                obj.ChebRoot = ChebRoot;
            end
            
            if nargin<6
                refine(obj,f,grid_opt);
            else
                obj.TreeGrid = obj.ChebRoot.leafGrids();
                
                obj.leafArray = obj.ChebRoot.collectLeaves();
            end
            
        end
        
        function refine(obj,f,grid_opt)
            
            if nargin<3
                grid_opt = false;
            end
            
            while ~obj.ChebRoot.is_refined
                
                Max = obj.ChebRoot.sample(f,grid_opt);
                
                if obj.ChebRoot.is_leaf
                    obj.ChebRoot = obj.ChebRoot.splitleaf(Max);
                else
                    obj.ChebRoot.PUsplit(Max);
                end
                
            end
            
            obj.TreeGrid = obj.ChebRoot.leafGrids();
            
            if ~obj.ChebRoot.is_leaf
                obj.ChebRoot.findIndex([]);
                obj.leafArray = obj.ChebRoot.collectLeaves();
            else
                obj.leafArray = {obj.ChebRoot};
            end
            
        end
        
        function ef = evalf(obj,X)
            ef = obj.ChebRoot.evalf(X);
        end
        
        function ef = evalfGrid(obj,X)
            ef = obj.ChebRoot.evalfGrid(X);
        end
        
        function ef = evalfZoneGrid(obj,X)
            if obj.ChebRoot.is_leaf
                ef = obj.ChebRoot.evalfGrid(X);
            else
                ef = obj.ChebRoot.evalfZoneGrid(X);
            end
        end
        
        function addTree = plus(obj,Tree2)
            
            add_f = @(x) obj.evalfGrid(x) + Tree2.evalfGrid(x);
            addTreeRoot = PUFun.add(obj.ChebRoot,Tree2.ChebRoot,add_f);
            
            addTree = PUFun(obj.domain,obj.deg_in,[],obj.tol,true,addTreeRoot);
            
        end
        
        function subTree = minus(obj,Tree2)
            
            sub_f = @(x) obj.evalfGrid(x) - Tree2.evalfGrid(x);
            subTreeRoot = PUFun.subtract(obj.ChebRoot,Tree2.ChebRoot,sub_f);
            
            subTree = PUFun(obj.domain,obj.deg_in,[],obj.tol,true,subTreeRoot);
            
        end
        
        function MultTree = mtimes(obj,Tree2)
            
            mult_f = @(x) obj.evalfGrid(x).*Tree2.evalfGrid(x);
            multTreeRoot = PUFun.multiply(obj.ChebRoot,Tree2.ChebRoot,mult_f);
            
            MultTree = PUFun(obj.domain,obj.deg_in,[],obj.tol,true,multTreeRoot);
            
        end
        
        function DivTree = mrdivide(obj,Tree2)
            
            div_f = @(x)obj.evalfGrid(x)./Tree2.evalfGrid(x);
            divTreeRoot = PUFun.divide(obj.ChebRoot,Tree2.ChebRoot,div_f);
            
            DivTree = PUFun(obj.domain,obj.deg_in,[],obj.tol,true,divTreeRoot);
            
        end
        
        function PowTree = mpower(obj,p)
            
            PowTreeRoot = PUFun.power(obj.ChebRoot,p);
            
            PowTree = PUFun(obj.domain,obj.deg_in,[],obj.tol,true,PowTreeRoot);
            
        end
        
        function ef = evalfTreeGrid(obj)
            
            for i=1:length(obj.TreeGrid)
                ef = obj.ChebRoot.evalfZoneGrid(obj.TreeGrid{i});
            end
            
        end
        
        function Coarsen(obj)
            obj.ChebRoot.Coarsen
        end
        
        
        % [PUF,DX,DXX]=evalf(obj,X)
        % This method evalutes the PUM approximation at X.
        %Input:
        %   x          : this constructs matrices given f|x.
        %Output:
        %   PU         : vector of approximation at X
        %   DX         : vector of derivative values at X
        %   DXX        : vector of second derivative values at X
        function int = sum(obj)
            
            int = 0;
            for i=1:length(obj.leafArray)
                
                if ~isempty(obj.leafArray{i}.zone)
                    X = cell(1,obj.ChebRoot.dim);
                    
                    for j=1:obj.ChebRoot.dim
                        [X{j},W{j}] = chebpts(obj.leafArray{i}.degs(j),obj.leafArray{i}.zone(j,:));
                    end
                    
                    vals = obj.leafArray{i}.evalfGrid(X);
                    
                    
                    if obj.ChebRoot.dim==2
                        int = int + W{2}*(W{1}*vals).';
                    else
                        int = int + W{2}*(W{1}*chebfun3.txm(vals,W{3},3))';
                    end
                    
                end
                
                
            end
        end
        
        function diff_Tree = diff(obj,diff_dim,order)
            diff_Tree = copy(obj);
            
            for i=1:length(obj.leafArray)
                obj.leafArray{i}.values = evalfDiffGrid(obj.leafArray{i},diff_dim,order);
            end
            
        end
        
        function ln = length(obj)
            ln = length(obj.ChebRoot);
        end
        
        function disp(obj)
            disp(obj.ChebRoot.toString());
        end
        
        function show(obj)
            show(obj.ChebRoot)
        end
        
        function plot(obj)
            domain = obj.ChebRoot.domain;
            x = linspace(domain(1,1),domain(1,2),100)';
            y = linspace(domain(2,1),domain(2,2),100)';
            ef = obj.evalfGrid({x y});
            [X,Y] = ndgrid(x,y);
            surf(X,Y,ef);
            xlabel('x');
            ylabel('y');
        end
        
    end
    
    methods(Access = protected)
        % Override copyElement method:
        function cpObj = copyElement(obj)
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            
            cpObj.ChebRoot = obj.ChebRoot.copy();
            
            if ~cpObj.ChebRoot.is_leaf
                cpObj.ChebRoot.findIndex([]);
                cpObj.leafArray = cpObj.ChebRoot.collectLeaves();
            else
                cpObj.leafArray = cpObj.ChebRoot;
            end
        end
    end
    
    methods (Static)
        %This method slices an array along a certain dimension
        function out = slice(A, ix, dim)
            subses = repmat({':'}, [1 ndims(A)]);
            subses{dim} = ix;
            out = A(subses{:});
        end
        
        function T_add = add(T_1,T_2,add_f)
            
            if T_1.is_leaf
                
                T_add = copy(T_2);
                
                leafArray = T_add.collectLeaves();
                
                for i=1:length(leafArray)
                    if isequal(T_1.domain,leafArray{i}.domain) && isequal(T_1.degs,leafArray{i}.degs)
                        leafArray{i}.values = leafArray{i}.values+T_1.values;
                    else
                        leafArray{i}.values = leafArray{i}.values+T_1.evalfGrid(leafArray{i}.leafGrids());
                    end
                end
                
            elseif T_2.is_leaf
                T_add = copy(T_1);
                
                leafArray = T_add.collectLeaves();
                
                for i=1:length(leafArray)
                    if isequal(T_2.domain,leafArray{i}.domain) && isequal(T_2.degs,leafArray{i}.degs)
                        leafArray{i}.values = leafArray{i}.values+T_2.values;
                    else
                        leafArray{i}.values = leafArray{i}.values+T_2.evalfGrid(leafArray{i}.leafGrids());
                    end
                end
                
            elseif T_1.splitting_dim == T_2.splitting_dim
                
                children{1} = PUFun.add(T_1.children{1},T_2.children{1},add_f);
                children{2} = PUFun.add(T_1.children{2},T_2.children{2},add_f);
                
                cheb_length_add = length(T_1.children{1})+length(T_1.children{2});
                domain_add = [T_1.children{1}.domain(:,1) T_1.children{2}.domain(:,2)];
                
                T_add = PUPatch(domain_add,T_1.zone,cheb_length_add,children,T_1.splitting_dim,T_1.index);
            else
                leafArray = T_1.collectLeaves();
                num_patch1 = length(leafArray);
                
                leafArray = T_2.collectLeaves();
                num_patch2 = length(leafArray);
                
                if num_patch1>num_patch2
                    T_add = copy(T_1);
                else
                    T_add = copy(T_2);
                end
                
                T_add.reset();
                T_add = T_add.refine(add_f,true);
            end
            
        end
        
        function T_sub = subtract(T_1,T_2,sub_f)
            
            if T_1.is_leaf
                
                T_sub = copy(T_2);
                
                if ~T_sub.is_leaf
                    T_sub.findIndex([]);
                    leafArray = T_sub.collectLeaves();
                else
                    leafArray = {T_sub};
                end
                
                for i=1:length(leafArray)
                    if isequal(T_1.domain,leafArray{i}.domain) && isequal(T_1.degs,leafArray{i}.degs)
                        leafArray{i}.values = T_1.values-leafArray{i}.values;
                    else
                        leafArray{i}.values = T_1.evalfGrid(leafArray{i}.leafGrids())-leafArray{i}.values;
                    end
                end
                
            elseif T_2.is_leaf
                
                T_sub = copy(T_1);
                
                leafArray = T_sub.collectLeaves();
                
                for i=1:length(leafArray)
                    if isequal(T_2.domain,leafArray{i}.domain) && isequal(T_2.degs,leafArray{i}.degs)
                        leafArray{i}.values = leafArray{i}.values-T_2.values;
                    else
                        leafArray{i}.values = leafArray{i}.values-T_2.evalfGrid(leafArray{i}.leafGrids());
                    end
                end
                
            elseif T_1.splitting_dim == T_2.splitting_dim
                
                children{1} = PUFun.subtract(T_1.children{1},T_2.children{1},sub_f);
                children{2} = PUFun.subtract(T_1.children{2},T_2.children{2},sub_f);
                
                cheb_length_add = length(T_1.children{1})+length(T_1.children{2});
                domain_add = [T_1.children{1}.domain(:,1) T_1.children{2}.domain(:,2)];
                
                T_sub = PUPatch(domain_add,T_1.zone,cheb_length_add,children,T_1.splitting_dim,T_1.index);
            else
                leafArray = T_1.collectLeaves();
                num_patch1 = length(leafArray);
                
                leafArray = T_2.collectLeaves();
                num_patch2 = length(leafArray);
                
                if num_patch1>num_patch2
                    T_sub = copy(T_1);
                else
                    T_sub = copy(T_2);
                end
                
                T_sub.reset();
                T_sub = T_sub.refine(sub_f,true);
            end
            
        end
        
        function T_mult = multiply(T_1,T_2,mult_f)
            
            if T_1.is_leaf
                
                T_mult = copy(T_2);
                
                leafArray = T_mult.collectLeaves();
                T_2Array = T_2.collectLeaves();
                
                T_mult.reset();
                
                for i=1:length(leafArray)
                    leafArray{i} = refine(leafArray{i},@(x)T_1.evalfGrid(x).*T_2Array{i}.evalfGrid(x),true);
                end
                
            elseif T_2.is_leaf
                T_mult = copy(T_1);
                
                
                T_Array = T_1.collectLeaves();
                leafArray = T_mult.collectLeaves();
                
                T_mult.reset();
                
                for i=1:length(leafArray)
                    leafArray{i} = refine(leafArray{i},@(x)T_2.evalfGrid(x).*T_Array{i}.evalfGrid(x),true);
                end
                
            elseif T_1.splitting_dim == T_2.splitting_dim
                
                children{1} = PUFun.multiply(T_1.children{1},T_2.children{1},mult_f);
                children{2} = PUFun.multiply(T_1.children{2},T_2.children{2},mult_f);
                
                cheb_length_add = length(T_1.children{1})+length(T_1.children{2});
                domain_add = [T_1.children{1}.domain(:,1) T_1.children{2}.domain(:,2)];
                
                T_mult = PUPatch(domain_add,T_1.zone,cheb_length_add,children,T_1.splitting_dim,T_1.index);
            else
                leafArray = T_1.collectLeaves();
                num_patch1 = length(leafArray);
                
                leafArray = T_2.collectLeaves();
                num_patch2 = length(leafArray);
                
                if num_patch1>num_patch2
                    T_mult = copy(T_1);
                else
                    T_mult = copy(T_2);
                end
                
                T_mult.reset();
                T_mult = T_mult.refine(mult_f,true);
            end
        end
        
        function T_divide = divide(T_1,T_2,divide_f)
            
            if T_1.is_leaf
                
                T_divide = copy(T_2);
                
                leafArray = T_divide.collectLeaves();
                T_2Array = T_2.collectLeaves();
                
                T_divide.reset();
                
                for i=1:length(leafArray)
                    leafArray{i} = refine(leafArray{i},@(x)T_1.evalfGrid(x)./T_2Array{i}.evalfGrid(x),true);
                end
                
            elseif T_2.is_leaf
                T_divide = copy(T_1);
                
                T_Array = T_1.collectLeaves();
                leafArray = T_divide.collectLeaves();
                
                T_divide.reset();
                
                for i=1:length(leafArray)
                    leafArray{i} = refine(leafArray{i},@(x)T_Array{i}.evalfGrid(x)./T_2.evalfGrid(x),true);
                end
                
            elseif T_1.splitting_dim == T_2.splitting_dim
                
                children{1} = PUFun.divide(T_1.children{1},T_2.children{1},divide_f);
                children{2} = PUFun.divide(T_1.children{2},T_2.children{2},divide_f);
                
                cheb_length_add = length(T_1.children{1})+length(T_1.children{2});
                domain_add = [T_1.children{1}.domain(:,1) T_1.children{2}.domain(:,2)];
                
                T_divide = PUPatch(domain_add,T_1.zone,cheb_length_add,children,T_1.splitting_dim,T_1.index);
            else
                
                leafArray = T_1.collectLeaves();
                num_patch1 = length(leafArray);
                
                leafArray = T_2.collectLeaves();
                num_patch2 = length(leafArray);
                
                if num_patch1>num_patch2
                    T_divide = copy(T_1);
                else
                    T_divide = copy(T_2);
                end
                
                T_divide.reset();
                T_divide = T_divide.refine(divide_f,true);
            end
        end
        
        function T_power = power(Tree,p)
            T_power = copy(Tree);
            
            T_power.reset();
            
            T_Array = Tree.collectLeaves();
            leafArray = T_power.collectLeaves();
            
            
            for i=1:length(T_Array)
                leafArray{i} = refine(leafArray{i},@(x)T_Array{i}.evalfGrid(x).^p,true);
            end
        end
    end
end

