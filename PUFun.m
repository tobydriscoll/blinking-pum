classdef PUFun < handle
    
    properties
        ChebRoot
        TreeGrid
        leafArray
    end
    
    methods
        
        function obj = PUFun(domain,deg_in,f,tol)
            
            if nargin < 4
                tol = 1e-12;
            end
                
            [dim,~] = size(domain);
            obj.ChebRoot = ChebPatch(domain,deg_in,ones(1,dim),tol);
            
            %Refine on f(x)
            
            while ~obj.ChebRoot.is_refined
                
                obj.ChebRoot.sample(f);
                
                if obj.ChebRoot.is_leaf
                    obj.ChebRoot = obj.ChebRoot.splitleaf();
                else
                    obj.ChebRoot.PUsplit();
                end
                
            end
            obj.TreeGrid = obj.ChebRoot.leafGrids();
            
            if ~obj.ChebRoot.is_leaf
                obj.ChebRoot.findIndex([]);
                obj.leafArray = obj.ChebRoot.collectLeaves({});
            else
                obj.leafArray = obj.ChebRoot;
            end
            
            
        end
        
        function ef = evalf(obj,X,dim,order)
            ef = obj.ChebRoot.evalf(X,dim,order);
        end
        
        function ef = evalfGrid(obj,X,dim,order)
            ef = obj.ChebRoot.evalfGrid(X,dim,order);
        end
        
        function ef = evalfTreeGrid(obj,dim,order)
            
            for i=1:length(obj.TreeGrid)
                ef = obj.ChebRoot.evalfGrid(obj.TreeGrid{i},dim,order);
            end
            
        end
        
        
        % [PUF,DX,DXX]=evalf(obj,X)
        % This method evalutes the PUM approximation at X.
        %Input:
        %   x          : this constructs matrices given f|x.
        %Output:
        %   PU         : vector of approximation at X
        %   DX         : vector of derivative values at X
        %   DXX        : vector of second derivative values at X
        function int = sum(obj,N)
            
            int = 0;
            for i=1:length(obj.leafArray)
                
                X = cell(1,obj.ChebRoot.dim);
                W = cell(1,obj.ChebRoot.dim);
                
                for j=1:obj.ChebRoot.dim
                    [X{j},W{j}] = chebpts(N,obj.leafArray{i}.domain(j,:));
                end
                
                WEIGHTSVALS = obj.ChebRoot.evalweights(obj.leafArray{i}.index,X,1,0);
                
                WEIGHTSG = cell(1,obj.ChebRoot.dim);
                [WEIGHTSG{:}] = ndgrid(WEIGHTSVALS{:});
                
                WEIGHTS = ones(size(WEIGHTSG{1}));
                
                for j=1:obj.ChebRoot.dim
                    WEIGHTS = WEIGHTS.*WEIGHTSG{j};
                end

                vals = WEIGHTS.*(obj.leafArray{i}.evalfGrid(X,1,0));
                
                if obj.ChebRoot.dim==3
                    vals = chebfun3t.unfold(vals,3);
                    vals = reshape(W{3}*vals,length(X{1}),length(X{2}));
                end
                
                int = int + W{2}*(W{1}*vals).';
                
            end
        end
        
        
        function ln = length(obj)
            ln = length(obj.ChebRoot);
        end
        
        function disp(obj)
            disp(obj.ChebRoot.toString());
        end
        
    end
end

