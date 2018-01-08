classdef PUFun < handle
    
    properties
        ChebRoot
        TreeGrid
        leafArray
        Errs
        Nums
    end
    
    methods
        
        function obj = PUFun(domain,deg_in,f,tol)
            
            if nargin < 4
                tol = 1e-12;
            end
                
            [dim,~] = size(domain);
            obj.ChebRoot = ChebPatch(domain,domain,domain,deg_in,true(1,dim),tol);
            
            %Refine on f(x)
            
            errs = [];
            
            nums = [];
%             x = linspace(domain(1,1),domain(1,2),100).';
%             y = linspace(domain(2,1),domain(2,2),100).';
%             
%             [X,Y] = ndgrid(x,y);
            while ~obj.ChebRoot.is_refined
                
                Max = obj.ChebRoot.sample(f);
                
%                 E = obj.ChebRoot.evalfGrid({x,y});
%                 
%                 errs = [errs norm(E(:)-f([X(:) Y(:)]),inf)];
%                 
%                 nums = [nums length(obj.ChebRoot)];
                
                if obj.ChebRoot.is_leaf
                    obj.ChebRoot = obj.ChebRoot.splitleaf(Max);
                else
                    obj.ChebRoot.PUsplit(Max);
                end
                
            end
            
            obj.Errs = errs;
            obj.Nums = nums;
            obj.TreeGrid = obj.ChebRoot.leafGrids();
            
            if ~obj.ChebRoot.is_leaf
                obj.ChebRoot.findIndex([]);
                obj.leafArray = obj.ChebRoot.collectLeaves({});
            else
                obj.leafArray = obj.ChebRoot;
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
        
        function ef = evalfTreeGrid(obj)
            
            for i=1:length(obj.TreeGrid)
                ef = obj.ChebRoot.evalfGrid(obj.TreeGrid{i});
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
end

