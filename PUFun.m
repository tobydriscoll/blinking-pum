classdef PUFun < handle & matlab.mixin.Copyable
    % This is the class for the PU approximation on squares and cubes.  
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
            obj.tol = tol;
            
            if nargin < 5
                grid_opt = false;
                obj.ChebRoot = ChebPatch(domain,domain,domain,deg_in,true(1,dim),obj.tol);
            elseif nargin < 6
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
        
        % refine(obj,f,grid_opt)
        % This method refines the tree by f(x).
        %Input:
        %   f          : the function to split on
        %   grid_opt   : boolean value indicating if
        %                function is evaluated for grids;
        %                must take cell array of grids
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
            
            obj.ChebRoot.clean();
            
        end
        
        % ef = evalf(obj,X)
        % This method evalutes the PUM approximation at X.
        %Input:
        %   x          : this constructs matrices given f|x.
        %Output:
        %   ef         : vector of approximation at X
        function ef = evalf(obj,X)
            ef = obj.ChebRoot.evalf(X);
        end
        
        % ef = evalf(obj,X)
        % This method evalutes the PUM approximation at at grid X.
        %Input:
        %   x          : this constructs matrices given f|x.
        %Output:
        %   ef         : grid of approximation at X
        function ef = evalfGrid(obj,X)
            ef = obj.ChebRoot.evalfGrid(X);
        end
        
        % addTree = plus(obj,Tree2)
        % This method adds obj and Tree2
        %Input:
        %   Tree2      : the other tree
        %Output:
        %   addTree    : new tree of the sum
        function addTree = plus(obj,Tree2)
            
            add_T = ChebPatch(obj.domain,obj.domain,obj.domain);
            
            add_T.deg_in = max(obj.deg_in,Tree2.deg_in);
            add_T.degs = add_T.standard_variables.chebpoints{add_T.deg_in};
            
            
            addTreeRoot = PUFun.fast_add(obj.ChebRoot,Tree2.ChebRoot,add_T,0);
            
            addTreeRoot.clean();
            
            addTree = PUFun(add_T.domain,add_T.deg_in,@(x)x,add_T.tol,false,addTreeRoot);
            
        end
        
        % subTree = minus(obj,Tree2)
        % This method subtracts obj and Tree2
        %Input:
        %   Tree2      : the other tree
        %Output:
        %   subTree    : new tree of the difference
        function subTree = minus(obj,Tree2)
            
            sub_T = ChebPatch(obj.domain,obj.domain,obj.domain);
            
            sub_T.deg_in = max(obj.deg_in,Tree2.deg_in);
            sub_T.degs = sub_T.standard_variables.chebpoints{sub_T.deg_in};
            
            
            subTreeRoot = PUFun.fast_add(obj.ChebRoot,Tree2.ChebRoot,sub_T,0);
            
            subTreeRoot.clean();
            
            subTree = PUFun(sub_T.domain,sub_T.deg_in,@(x)x,sub_T.tol,false,subTreeRoot);
        end
        
        % MultTree = mtimes(obj,Tree2)
        % This method multiplies obj and Tree2
        %Input:
        %   Tree2      : the other tree
        %Output:
        %   MultTree   : new tree of the product
        function MultTree = mtimes(obj,Tree2)
            
            mult_T = ChebPatch(obj.domain,obj.domain,obj.domain);
            
            mult_T.deg_in = max(obj.deg_in,Tree2.deg_in);
            mult_T.degs = mult_T.standard_variables.chebpoints{mult_T.deg_in};
            
            
            multTreeRoot = PUFun.fast_multiply(obj.ChebRoot,Tree2.ChebRoot,mult_T,0);
            
            multTreeRoot.clean();
            
            MultTree = PUFun(mult_T.domain,mult_T.deg_in,@(x)x,mult_T.tol,false,multTreeRoot);
            
        end
        
        % DivTree = mrdivide(obj,Tree2)
        % This method divides obj and Tree2
        %Input:
        %   Tree2      : the other tree
        %Output:
        %   DivTree    : new tree of the quotient
        function DivTree = mrdivide(obj,Tree2)
            
            div_T = ChebPatch(obj.domain,obj.domain,obj.domain);
            
            div_T.deg_in = max(obj.deg_in,Tree2.deg_in);
            div_T.degs = div_T.standard_variables.chebpoints{div_T.deg_in};
            
            
            DivTreeRoot = PUFun.fast_multiply(obj.ChebRoot,Tree2.ChebRoot,div_T,0);
            
            DivTreeRoot.clean();
            
            DivTree = PUFun(div_T.domain,div_T.deg_in,@(x)x,div_T.tol,false,DivTreeRoot);
            
            
        end
        
        % PowTree = mpower(obj,p)
        % This method computes obj to the power p
        %Input:
        %   p          : the power to be used
        %Output:
        %   PowTree    : the new tree of the power
        function PowTree = mpower(obj,p)
            
            PowTreeRoot = PUFun.power(obj.ChebRoot,p);
            
            PowTree = PUFun(obj.domain,obj.deg_in,[],obj.tol,true,PowTreeRoot);
            
        end
        
        % Coarsen(obj)
        % This method Coarsens each of the patches
        function Coarsen(obj)
            obj.ChebRoot.Coarsen
        end
        
        % Refines(obj)
        % This method Refines each of the patches
        function Refine(obj)
            obj.ChebRoot.Coarsen
        end
        
        
        % int = sum(obj)
        % This method computes the integral of obj
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
        
        % diff_Tree = diff(obj,diff_dim,order)
        % This method computes thd approximation of the derivative
        %Input:
        %Output:
        %   diff_Tree  : PU approximation of derivative
        function diff_Tree = diff(obj,diff_dim,order)
            diff_Tree = copy(obj);
            
            for i=1:length(obj.leafArray)
                diff_Tree.leafArray{i}.values = evalfDiffGrid(obj.leafArray{i},diff_dim,order);
            end
            
        end
        
        % div_Tree = div(obj)
        % This method computes thd approximation of the derivative
        %Input:
        %Output:
        %   div_Tree  : PU approximation of divergence
        function div_Tree = div(obj)
            div_Tree = copy(obj);
            
            for i=1:length(obj.leafArray)
                if obj.ChebRoot.dim==1
                    div_Tree.leafArray{i}.values = evalfDiffGrid(obj.leafArray{i},1,1);
                elseif obj.ChebRoot.dim==2
                    div_Tree.leafArray{i}.values = evalfDiffGrid(obj.leafArray{i},1,1)+evalfDiffGrid(obj.leafArray{i},2,1);
                else
                    div_Tree.leafArray{i}.values = evalfDiffGrid(obj.leafArray{i},1,1)+evalfDiffGrid(obj.leafArray{i},2,1)+evalfDiffGrid(obj.leafArray{i},3,1);
                end
            end
            
        end
        
        % lap_Tree = lap(obj)
        % This method computes thd approximation of the Laplacian
        %Input:
        %Output:
        %   lap_Tree  : PU approximation of Laplacian
        function lap_Tree = lap(obj)
            lap_Tree = copy(obj);
            
            for i=1:length(obj.leafArray)
                if obj.ChebRoot.dim==1
                    lap_Tree.leafArray{i}.values = evalfDiffGrid(obj.leafArray{i},1,2);
                elseif obj.ChebRoot.dim==2
                    lap_Tree.leafArray{i}.values = evalfDiffGrid(obj.leafArray{i},1,2)+evalfDiffGrid(obj.leafArray{i},2,2);
                else
                    lap_Tree.leafArray{i}.values = evalfDiffGrid(obj.leafArray{i},1,2)+evalfDiffGrid(obj.leafArray{i},2,2)+evalfDiffGrid(obj.leafArray{i},3,2);
                end
            end
            
        end
        
        % ln = length(obj)
        % This returns the length of the tree
        function ln = length(obj)
            ln = length(obj.ChebRoot);
        end
        
        % disp(obj)
        % Returns string of object
        function disp(obj)
            disp(obj.ChebRoot.toString());
        end
        
        % disp(obj)
        % Returns color plot of patches
        function show(obj)
            show(obj.ChebRoot)
        end
        
        % plot(obj)
        % Returns a plot of the function in 2D
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
            
            cpObj.leafArray = cpObj.ChebRoot.collectLeaves();
        end
    end
    
    methods (Static)
        %This method slices an array along a certain dimension
        function out = slice(A, ix, dim)
            subses = repmat({':'}, [1 ndims(A)]);
            subses{dim} = ix;
            out = A(subses{:});
        end
        
        
        % T_add = add(T_1,T_2,add_f)
        % This method adds T_1 and T_2
        %Input:
        %   T_1,T_2    : trees to be added together
        %   add_f      : function for adding T_1 and T_2 using PU
        %Output:
        %   T_add      : PU approximation of the sum
        function T_add = fast_add(T_1,T_2,T_add,split_dim)
            
%             if ~any(T_1.zone(:,1)<=T_add.zone(:,1) & T_1.zone(:,2)>=T_add.zone(:,2)) | ~any(T_1.zone(:,1)<=T_add.zone(:,1) & T_1.zone(:,1)>=T_add.zone(:,1))
%                 error('the domains dont work');
%             end
            
            if T_1.is_leaf && T_2.is_leaf
                T_add.deg_in = max(T_1.deg_in,T_2.deg_in);
                T_add.degs = max(T_1.degs,T_2.degs);
                T_add.values = T_1.evalfGrid(T_add.leafGrids())+T_2.evalfGrid(T_add.leafGrids());
                %T_add = T_add.refine(@(x)T_1.evalfGrid(x)+T_2.evalfGrid(x),true);
            elseif T_1.is_leaf && ~T_2.is_leaf
                T_add = T_add.split(T_2.splitting_dim);
                T_add.children{1} = PUFun.fast_add(T_1,T_2.children{1},T_add.children{1},T_add.splitting_dim);
                T_add.children{2} = PUFun.fast_add(T_1,T_2.children{2},T_add.children{2},T_add.splitting_dim);
            elseif ~T_1.is_leaf && T_2.is_leaf
                T_add = T_add.split(T_1.splitting_dim);
                T_add.children{1} = PUFun.fast_add(T_1.children{1},T_2,T_add.children{1},T_add.splitting_dim);
                T_add.children{2} = PUFun.fast_add(T_1.children{2},T_2,T_add.children{2},T_add.splitting_dim);
            elseif T_1.splitting_dim == T_2.splitting_dim
                T_add = T_add.split(T_1.splitting_dim);
                T_add.children{1} = PUFun.fast_add(T_1.children{1},T_2.children{1},T_add.children{1},T_add.splitting_dim);
                T_add.children{2} = PUFun.fast_add(T_1.children{2},T_2.children{2},T_add.children{2},T_add.splitting_dim);
            else
                
                if split_dim==0
                    next_split=1;
                else
                    next_split = mod(split_dim,T_1.dim)+1;
                end
                
                while true
                    if T_1.splitting_dim == next_split
                        First_T = T_1;
                        Last_T = T_2;
                        break;
                    elseif T_2.splitting_dim == next_split
                        First_T = T_2;
                        Last_T = T_1;
                        break;
                    else
                        next_split = mod(next_split,T_1.dim)+1;
                    end
                end
                
                
                if ~Last_T.split_flag(First_T.splitting_dim)
                    T_add = T_add.split(First_T.splitting_dim);
                    T_add.children{1} = PUFun.fast_add(First_T.children{1},Last_T,T_add.children{1},T_add.splitting_dim);
                    T_add.children{2} = PUFun.fast_add(First_T.children{2},Last_T,T_add.children{2},T_add.splitting_dim);
                else
                    T_add = T_add.refine(@(x)T_1.evalfGrid(x)+T_2.evalfGrid(x),true);
                end
                
                
            end
        end
        
                % T_add = add(T_1,T_2,add_f)
        % This method adds T_1 and T_2
        %Input:
        %   T_1,T_2    : trees to be added together
        %   add_f      : function for adding T_1 and T_2 using PU
        %Output:
        %   T_add      : PU approximation of the sum
        function T_sub = fast_subtract(T_1,T_2,T_sub,split_dim)
            
            
            if T_1.is_leaf && T_2.is_leaf
                T_sub.deg_in = max(T_1.deg_in,T_2.deg_in);
                T_sub.degs = max(T_1.degs,T_2.degs);
                T_sub.values = T_1.evalfGrid(T_sub.leafGrids())-T_2.evalfGrid(T_sub.leafGrids());
            elseif T_1.is_leaf && ~T_2.is_leaf
                T_sub = T_sub.split(T_2.splitting_dim);
                T_sub.children{1} = PUFun.fast_subtract(T_1,T_2.children{1},T_sub.children{1},T_sub.splitting_dim);
                T_sub.children{2} = PUFun.fast_subtract(T_1,T_2.children{2},T_sub.children{2},T_sub.splitting_dim);
            elseif ~T_1.is_leaf && T_2.is_leaf
                T_sub = T_sub.split(T_1.splitting_dim);
                T_sub.children{1} = PUFun.fast_subtract(T_1.children{1},T_2,T_sub.children{1},T_sub.splitting_dim);
                T_sub.children{2} = PUFun.fast_subtract(T_1.children{2},T_2,T_sub.children{2},T_sub.splitting_dim);
            elseif T_1.splitting_dim == T_2.splitting_dim
                T_sub = T_sub.split(T_1.splitting_dim);
                T_sub.children{1} = PUFun.fast_subtract(T_1.children{1},T_2.children{1},T_sub.children{1},T_sub.splitting_dim);
                T_sub.children{2} = PUFun.fast_subtract(T_1.children{2},T_2.children{2},T_sub.children{2},T_sub.splitting_dim);
            else
                
                if split_dim==0
                    next_split=1;
                else
                    next_split = mod(split_dim,T_1.dim)+1;
                end
                
                T_1_is_first = true;
                
                while true
                    if T_1.splitting_dim == next_split
                        T_1_is_first = true;
                        break;
                    elseif T_2.splitting_dim == next_split
                        T_1_is_first = false;
                        break;
                    else
                        next_split = mod(next_split,T_1.dim)+1;
                    end
                end
                
                
                if T_1_is_first && ~T_2.split_flag(T_1.splitting_dim)
                    
                    T_sub = T_sub.split(T_1.splitting_dim);
                    T_sub.children{1} = PUFun.fast_subtract(T_1.children{1},T_2,T_sub.children{1},T_sub.splitting_dim);
                    T_sub.children{2} = PUFun.fast_subtract(T_1.children{2},T_2,T_sub.children{2},T_sub.splitting_dim);
                    
                elseif T_2_is_first && ~T_1.split_flag(T_2.splitting_dim)
                    
                    T_sub = T_sub.split(T_2.splitting_dim);
                    T_sub.children{1} = PUFun.fast_subtract(T_1,T_2.children{1},T_sub.children{1},T_sub.splitting_dim);
                    T_sub.children{2} = PUFun.fast_subtract(T_1,T_2.children{2},T_sub.children{2},T_sub.splitting_dim);
                    
                else
                    T_sub = T_sub.refine(@(x)T_1.evalfGrid(x)-T_2.evalfGrid(x),true);
                end
                
                
            end
        end
        
        
                % T_add = add(T_1,T_2,add_f)
        % This method adds T_1 and T_2
        %Input:
        %   T_1,T_2    : trees to be added together
        %   add_f      : function for adding T_1 and T_2 using PU
        %Output:
        %   T_add      : PU approximation of the sum
        function T_mult = fast_multiply(T_1,T_2,T_mult,split_dim)
            
            
            if T_1.is_leaf && T_2.is_leaf
                T_mult = T_mult.refine(@(x)T_1.evalfGrid(x).*T_2.evalfGrid(x),true);
            elseif T_1.is_leaf && ~T_2.is_leaf
                T_mult = T_mult.split(T_2.splitting_dim);
                T_mult.children{1} = PUFun.fast_multiply(T_1,T_2.children{1},T_mult.children{1},T_mult.splitting_dim);
                T_mult.children{2} = PUFun.fast_multiply(T_1,T_2.children{2},T_mult.children{2},T_mult.splitting_dim);
            elseif ~T_1.is_leaf && T_2.is_leaf
                T_mult = T_mult.split(T_1.splitting_dim);
                T_mult.children{1} = PUFun.fast_multiply(T_1.children{1},T_2,T_mult.children{1},T_mult.splitting_dim);
                T_mult.children{2} = PUFun.fast_multiply(T_1.children{2},T_2,T_mult.children{2},T_mult.splitting_dim);
            elseif T_1.splitting_dim == T_2.splitting_dim
                T_mult = T_mult.split(T_1.splitting_dim);
                T_mult.children{1} = PUFun.fast_multiply(T_1.children{1},T_2.children{1},T_mult.children{1},T_mult.splitting_dim);
                T_mult.children{2} = PUFun.fast_multiply(T_1.children{2},T_2.children{2},T_mult.children{2},T_mult.splitting_dim);
            else
                
                if split_dim==0
                    next_split=1;
                else
                    next_split = mod(split_dim,T_1.dim)+1;
                end
                
                while true
                    if T_1.splitting_dim == next_split
                        First_T = T_1;
                        Last_T = T_2;
                        break;
                    elseif T_2.splitting_dim == next_split
                        First_T = T_2;
                        Last_T = T_1;
                        break;
                    else
                        next_split = mod(next_split,T_1.dim)+1;
                    end
                end
                
                
                if ~Last_T.split_flag(First_T.splitting_dim)
                    T_mult = T_mult.split(First_T.splitting_dim);
                    T_mult.children{1} = PUFun.fast_multiply(First_T.children{1},Last_T,T_mult.children{1},T_mult.splitting_dim);
                    T_mult.children{2} = PUFun.fast_multiply(First_T.children{2},Last_T,T_mult.children{2},T_mult.splitting_dim);
                else
                    T_mult = T_mult.refine(@(x)T_1.evalfGrid(x).*T_2.evalfGrid(x),true);
                end
                
                
            end
        end
        
                
                        % T_add = add(T_1,T_2,add_f)
        % This method adds T_1 and T_2
        %Input:
        %   T_1,T_2    : trees to be added together
        %   add_f      : function for adding T_1 and T_2 using PU
        %Output:
        %   T_add      : PU approximation of the sum
        function T_div = fast_divide(T_1,T_2,T_div,split_dim)
            
            
            if T_1.is_leaf && T_2.is_leaf
                T_div = T_div.refine(@(x)T_1.evalfGrid(x)./T_2.evalfGrid(x),true);
            elseif T_1.is_leaf && ~T_2.is_leaf
                T_div = T_div.split(T_2.splitting_dim);
                T_div.children{1} = PUFun.fast_divide(T_1,T_2.children{1},T_div.children{1},T_div.splitting_dim);
                T_div.children{2} = PUFun.fast_divide(T_1,T_2.children{2},T_div.children{2},T_div.splitting_dim);
            elseif ~T_1.is_leaf && T_2.is_leaf
                T_div = T_div.split(T_1.splitting_dim);
                T_div.children{1} = PUFun.fast_divide(T_1.children{1},T_2,T_div.children{1},T_div.splitting_dim);
                T_div.children{2} = PUFun.fast_divide(T_1.children{2},T_2,T_div.children{2},T_div.splitting_dim);
            elseif T_1.splitting_dim == T_2.splitting_dim
                T_div = T_div.split(T_1.splitting_dim);
                T_div.children{1} = PUFun.fast_divide(T_1.children{1},T_2.children{1},T_div.children{1},T_div.splitting_dim);
                T_div.children{2} = PUFun.fast_divide(T_1.children{2},T_2.children{2},T_div.children{2},T_div.splitting_dim);
            else
                
                if split_dim==0
                    next_split=1;
                else
                    next_split = mod(split_dim,T_1.dim)+1;
                end
                
                T_1_is_first = true;
                
                while true
                    if T_1.splitting_dim == next_split
                        T_1_is_first = true;
                        break;
                    elseif T_2.splitting_dim == next_split
                        T_1_is_first = false;
                        break;
                    else
                        next_split = mod(next_split,T_1.dim)+1;
                    end
                end
                
                
                if T_1_is_first && ~T_2.split_flag(T_1.splitting_dim)
                    
                    T_div = T_div.split(T_1.splitting_dim);
                    T_div.children{1} = PUFun.fast_divide(T_1.children{1},T_2,T_div.children{1},T_div.splitting_dim);
                    T_div.children{2} = PUFun.fast_divide(T_1.children{2},T_2,T_div.children{2},T_div.splitting_dim);
                    
                elseif T_2_is_first && ~T_1.split_flag(T_2.splitting_dim)
                    
                    T_div = T_div.split(T_2.splitting_dim);
                    T_div.children{1} = PUFun.fast_divide(T_1,T_2.children{1},T_div.children{1},T_div.splitting_dim);
                    T_div.children{2} = PUFun.fast_divide(T_1,T_2.children{2},T_div.children{2},T_div.splitting_dim);
                    
                else
                    T_div = T_div.refine(@(x)T_1.evalfGrid(x)./T_2.evalfGrid(x),true);
                end
                
                
            end
        end
        
        % T_power = power(Tree,p)
        % This method computes Tree^p
        %Input:
        %   Tree       : trees to be powered
        %   p          : power to use
        %Output:
        %   T_power    : PU approximation of the power
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


%         % T_add = add(T_1,T_2,add_f)
%         % This method adds T_1 and T_2
%         %Input:
%         %   T_1,T_2    : trees to be added together
%         %   add_f      : function for adding T_1 and T_2 using PU
%         %Output:
%         %   T_add      : PU approximation of the sum
%         function T_add = add(T_1,T_2,add_f)
%             
%             if T_1.is_leaf
%                 
%                 T_add = copy(T_2);
%                 
%                 leafArray = T_add.collectLeaves();
%                 
%                 for i=1:length(leafArray)
%                     if isequal(T_1.domain,leafArray{i}.domain) && isequal(T_1.degs,leafArray{i}.degs)
%                         leafArray{i}.values = leafArray{i}.values+T_1.values;
%                     else
%                         leafArray{i}.values = leafArray{i}.values+T_1.evalfGrid(leafArray{i}.leafGrids());
%                     end
%                 end
%                 
%             elseif T_2.is_leaf
%                 T_add = copy(T_1);
%                 
%                 leafArray = T_add.collectLeaves();
%                 
%                 for i=1:length(leafArray)
%                     if isequal(T_2.domain,leafArray{i}.domain) && isequal(T_2.degs,leafArray{i}.degs)
%                         leafArray{i}.values = leafArray{i}.values+T_2.values;
%                     else
%                         leafArray{i}.values = leafArray{i}.values+T_2.evalfGrid(leafArray{i}.leafGrids());
%                     end
%                 end
%                 
%             elseif T_1.splitting_dim == T_2.splitting_dim
%                 
%                 children{1} = PUFun.add(T_1.children{1},T_2.children{1},add_f);
%                 children{2} = PUFun.add(T_1.children{2},T_2.children{2},add_f);
%                 
%                 cheb_length_add = length(T_1.children{1})+length(T_1.children{2});
%                 domain_add = [T_1.children{1}.domain(:,1) T_1.children{2}.domain(:,2)];
%                 
%                 T_add = PUPatch(domain_add,T_1.zone,cheb_length_add,children,T_1.splitting_dim,T_1.index);
%             else
%                 leafArray = T_1.collectLeaves();
%                 num_patch1 = length(leafArray);
%                 
%                 leafArray = T_2.collectLeaves();
%                 num_patch2 = length(leafArray);
%                 
%                 if num_patch1>num_patch2
%                     T_add = copy(T_1);
%                 else
%                     T_add = copy(T_2);
%                 end
%                 
%                 T_add.reset();
%                 T_add = T_add.refine(add_f,true);
%             end
%             
%         end
        
%         % T_sub = subtract(T_1,T_2,sub_f)
%         % This method subtracts T_1 and T_2
%         %Input:
%         %   T_1,T_2    : trees to be subtracted together
%         %   sub_f      : function for subtracting T_1 and T_2 using PU
%         %Output:
%         %   T_sub      : PU approximation of the sum
%         function T_sub = subtract(T_1,T_2,sub_f)
%             
%             if T_1.is_leaf
%                 
%                 T_sub = copy(T_2);
%                 
%                 if ~T_sub.is_leaf
%                     T_sub.findIndex([]);
%                     leafArray = T_sub.collectLeaves();
%                 else
%                     leafArray = {T_sub};
%                 end
%                 
%                 for i=1:length(leafArray)
%                     if isequal(T_1.domain,leafArray{i}.domain) && isequal(T_1.degs,leafArray{i}.degs)
%                         leafArray{i}.values = T_1.values-leafArray{i}.values;
%                     else
%                         leafArray{i}.values = T_1.evalfGrid(leafArray{i}.leafGrids())-leafArray{i}.values;
%                     end
%                 end
%                 
%             elseif T_2.is_leaf
%                 
%                 T_sub = copy(T_1);
%                 
%                 leafArray = T_sub.collectLeaves();
%                 
%                 for i=1:length(leafArray)
%                     if isequal(T_2.domain,leafArray{i}.domain) && isequal(T_2.degs,leafArray{i}.degs)
%                         leafArray{i}.values = leafArray{i}.values-T_2.values;
%                     else
%                         leafArray{i}.values = leafArray{i}.values-T_2.evalfGrid(leafArray{i}.leafGrids());
%                     end
%                 end
%                 
%             elseif T_1.splitting_dim == T_2.splitting_dim
%                 
%                 children{1} = PUFun.subtract(T_1.children{1},T_2.children{1},sub_f);
%                 children{2} = PUFun.subtract(T_1.children{2},T_2.children{2},sub_f);
%                 
%                 cheb_length_add = length(T_1.children{1})+length(T_1.children{2});
%                 domain_add = [T_1.children{1}.domain(:,1) T_1.children{2}.domain(:,2)];
%                 
%                 T_sub = PUPatch(domain_add,T_1.zone,cheb_length_add,children,T_1.splitting_dim,T_1.index);
%             else
%                 leafArray = T_1.collectLeaves();
%                 num_patch1 = length(leafArray);
%                 
%                 leafArray = T_2.collectLeaves();
%                 num_patch2 = length(leafArray);
%                 
%                 if num_patch1>num_patch2
%                     T_sub = copy(T_1);
%                 else
%                     T_sub = copy(T_2);
%                 end
%                 
%                 T_sub.reset();
%                 T_sub = T_sub.refine(sub_f,true);
%             end
%             
%         end
        
        % T_mult = multiply(T_1,T_2,mult_f)
        % This method multiplies T_1 and T_2
        %Input:
        %   T_1,T_2    : trees to be multiplied together
        %   mult_f     : function for multiplying T_1 and T_2 using PU
        %Output:
        %   T_mult     : PU approximation of the product
%         function T_mult = multiply(T_1,T_2,mult_f)
%             
%             if T_1.is_leaf
%                 
%                 T_mult = copy(T_2);
%                 
%                 leafArray = T_mult.collectLeaves();
%                 T_2Array = T_2.collectLeaves();
%                 
%                 T_mult.reset();
%                 
%                 for i=1:length(leafArray)
%                     leafArray{i} = refine(leafArray{i},@(x)T_1.evalfGrid(x).*T_2Array{i}.evalfGrid(x),true);
%                 end
%                 
%             elseif T_2.is_leaf
%                 T_mult = copy(T_1);
%                 
%                 
%                 T_Array = T_1.collectLeaves();
%                 leafArray = T_mult.collectLeaves();
%                 
%                 T_mult.reset();
%                 
%                 for i=1:length(leafArray)
%                     leafArray{i} = refine(leafArray{i},@(x)T_2.evalfGrid(x).*T_Array{i}.evalfGrid(x),true);
%                 end
%                 
%             elseif T_1.splitting_dim == T_2.splitting_dim
%                 
%                 children{1} = PUFun.multiply(T_1.children{1},T_2.children{1},mult_f);
%                 children{2} = PUFun.multiply(T_1.children{2},T_2.children{2},mult_f);
%                 
%                 cheb_length_add = length(T_1.children{1})+length(T_1.children{2});
%                 domain_add = [T_1.children{1}.domain(:,1) T_1.children{2}.domain(:,2)];
%                 
%                 T_mult = PUPatch(domain_add,T_1.zone,cheb_length_add,children,T_1.splitting_dim,T_1.index);
%             else
%                 leafArray = T_1.collectLeaves();
%                 num_patch1 = length(leafArray);
%                 
%                 leafArray = T_2.collectLeaves();
%                 num_patch2 = length(leafArray);
%                 
%                 if num_patch1>num_patch2
%                     T_mult = copy(T_1);
%                 else
%                     T_mult = copy(T_2);
%                 end
%                 
%                 T_mult.reset();
%                 T_mult = T_mult.refine(mult_f,true);
%             end
%         end

        % T_divide = divide(T_1,T_2,divide_f)
        % This method divides T_1 and T_2
        %Input:
        %   T_1,T_2    : trees to be divided together
        %   mult_f     : function for dividing T_1 and T_2 using PU
        %Output:
        %   T_divide   : PU approximation of the quotient
%         function T_divide = divide(T_1,T_2,divide_f)
%             
%             if T_1.is_leaf
%                 
%                 T_divide = copy(T_2);
%                 
%                 leafArray = T_divide.collectLeaves();
%                 T_2Array = T_2.collectLeaves();
%                 
%                 T_divide.reset();
%                 
%                 for i=1:length(leafArray)
%                     leafArray{i} = refine(leafArray{i},@(x)T_1.evalfGrid(x)./T_2Array{i}.evalfGrid(x),true);
%                 end
%                 
%             elseif T_2.is_leaf
%                 T_divide = copy(T_1);
%                 
%                 T_Array = T_1.collectLeaves();
%                 leafArray = T_divide.collectLeaves();
%                 
%                 T_divide.reset();
%                 
%                 for i=1:length(leafArray)
%                     leafArray{i} = refine(leafArray{i},@(x)T_Array{i}.evalfGrid(x)./T_2.evalfGrid(x),true);
%                 end
%                 
%             elseif T_1.splitting_dim == T_2.splitting_dim
%                 
%                 children{1} = PUFun.divide(T_1.children{1},T_2.children{1},divide_f);
%                 children{2} = PUFun.divide(T_1.children{2},T_2.children{2},divide_f);
%                 
%                 cheb_length_add = length(T_1.children{1})+length(T_1.children{2});
%                 domain_add = [T_1.children{1}.domain(:,1) T_1.children{2}.domain(:,2)];
%                 
%                 T_divide = PUPatch(domain_add,T_1.zone,cheb_length_add,children,T_1.splitting_dim,T_1.index);
%             else
%                 
%                 leafArray = T_1.collectLeaves();
%                 num_patch1 = length(leafArray);
%                 
%                 leafArray = T_2.collectLeaves();
%                 num_patch2 = length(leafArray);
%                 
%                 if num_patch1>num_patch2
%                     T_divide = copy(T_1);
%                 else
%                     T_divide = copy(T_2);
%                 end
%                 
%                 T_divide.reset();
%                 T_divide = T_divide.refine(divide_f,true);
%             end
%         end

