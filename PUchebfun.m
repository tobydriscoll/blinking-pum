classdef PUchebfun < handle & matlab.mixin.Copyable
% PUFun   PUFun class for representing n-d functions on lines, rectangles
% and cuboids using partition of unitities.
% 
% This class represents smooth multivariate functions on hypercubes up to 
% dimension 3 with a partition of unity approximation. This class 
% automatically finds a set of overlapping domains that are adapted to the 
% features of the function, and the blends locally defined Chebyshev 
% approximations on the domains with a partition of unity.
% 
% PUFun(f) constructs a partition of unity approximation representing f on 
% the domain [-1 1]^k, where k is the number of variables of f 
% (i.e. the dimension of the domain of f). Functions must be vectorized.
% 
% PUFun(f,[a_1 b_1;a_2 b_2;...;a_d b_d]) constructs a partition of unity 
% approximation representing f on the domain [a_1 b_1;a_2 b_2;...;a_d b_d]. 
% Functions must be vectorized. 
% 
% PUFun(f,[a_1 b_1;a_2 b_2;...;a_d b_d],varargin) constructs a partition of unity approximation 
% representing f, based on the options passed into with varargin; that is 
% PUFun(f,'perf1',perf1,'pref2',pref2,..) is called. This preferences that 
% can be set are:
% 
% *The domain of the function: 'domain' , [a_1 b_1;a_2 b_2;...;a_d b_d]
% 
% *The degree indices from the standard degrees in each dimension : 
% 'degreeIndex', [ind_1,ind_2, ... ind_d]. 
% 
% % *The tolerance used for adaptation : 
% 'tol', tol. 
%
% Here the degrees can be chosen from the set [3 5 9 17 33 65 129].  
% So if 'degreeIndex', [5 5 5], the max degree of any approximate will be 
% 33 in each direction. 

    properties
        ChebRoot
        TreeGrid
        leafArray
        Errs
        Nums
        deg_in
        tol
        domain = [];
        grid_opt = false;
    end
    
    methods
        
        %function obj = PUFun(domain,deg_in,f,tol,grid_opt)
        function obj = PUchebfun(varargin)
            
                
                if length(varargin)==1
                    f = varargin{1};
                    obj.domain = repmat([-1 1],nargin(f),1);
                    cheb_struct.domain = obj.domain;
                    obj.ChebRoot = ChebPatch(cheb_struct);
                    
                elseif length(varargin)==2
                    f = varargin{1};
                    obj.domain = varargin{2};
                    cheb_struct.domain = obj.domain;
                    obj.ChebRoot = ChebPatch(cheb_struct);
                    
                else
                    f = varargin{1};
                    obj.domain = varargin{2};
                    cheb_struct.domain = obj.domain;
                    varargin(1:2) = [];
                    args = varargin;
                    while ( ~isempty(args) )
                        if strcmpi(args{1}, 'gridOption')
                            obj.grid_opt = args{2};
                        elseif strcmpi(args{1}, 'degreeIndex')
                            obj.deg_in = args{2};
                            cheb_struct.deg_in = args{2};
                        elseif strcmpi(args{1}, 'tol')
                            obj.tol = args{2};
                            cheb_struct.tol = args{2};
                        elseif strcmpi(args{1}, 'CourseDegreeIndex')
                            cheb_struct.cdeg_in = args{2};
                        else
                            error(strcat(args{1},' is not a valid parameter.'));
                        end
                        args(1:2) = [];
                    end
                end
                
                obj.cheb_deg_in = obj.ChebRoot.cheb_deg_in;
                obj.deg_in = obj.ChebRoot.deg_in;
                obj.tol = obj.ChebRoot.tol;
                
                obj.ChebRoot = ChebPatch(cheb_struct);
                refine(obj,f);
                
    end
            
        
        % refine(obj,f,grid_opt)
        % This method refines the tree by f(x).
        %Input:
        %   f          : the function to split on
        %   grid_opt   : boolean value indicating if
        %                function is evaluated for grids;
        %                must take cell array of grids
        function refine(obj,f)
            
            
            while ~obj.ChebRoot.is_refined
                
                Max = obj.ChebRoot.sample(f,obj.grid_opt);
                
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

