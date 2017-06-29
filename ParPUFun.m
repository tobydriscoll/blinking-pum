classdef ParPUFun < handle
    
    properties
        ChebRoot
        TreeGrid
        leafArray
    end
    
    methods
        
        function obj = ParPUFun(domain,deg_in,f,tol)
            
            if nargin < 4
                tol = 1e-12;
            end
                
            [dim,~] = size(domain);
            obj.ChebRoot = ChebPatch(domain,deg_in,ones(1,dim),tol);
            
           NewGen = {obj.ChebRoot};
            %Refine on f(x)
            
            while ~obj.ChebRoot.is_refined
                
                
                num_patches = ones(length(NewGen),1);
                
                %Determine if the leaves are refined in parallel
                parfor i=1:length(NewGen)
                    if ~NewGen{i}.is_refined
                        NewGen{i}.sample(f);
                        NewGen{i} = NewGen{i}.splitleaf();
                        if ~NewGen{i}.is_leaf
                            num_patches(i)=2;
                        end
                    end
                end
                
                index = 1;
                
                % The patches in PreviousGen are linked to ChebRoot. We
                % need to loop through and link the new children
                % to the tree rooted at ChebRoot (hence the names
                % new and previous generation).
                if length(NewGen)==1
                    obj.ChebRoot = NewGen{1};
                else
                    for i=1:length(PreviousGen)
                        if ~PreviousGen{i}.is_leaf
                            PreviousGen{i}.children{1} = NewGen{index};
                            PreviousGen{i}.children{2} = NewGen{index+1};
                            index = index + 2;
                        else
                            index = index + 1;
                        end
                    end
                end
                
                %Construct new array, update tree.
                newleaves = cell(sum(num_patches),1);
                index = 1;
                for i=1:length(NewGen)
                    if NewGen{i}.is_leaf
                        newleaves{index} = NewGen{i};
                        index = index + 1;
                    else
                        newleaves{index} = NewGen{i}.children{1};
                        newleaves{index+1} = NewGen{i}.children{2};
                        index = index + 2;
                    end
                end

                obj.ChebRoot.is_refined = length(NewGen)==length(newleaves);
                
                PreviousGen = NewGen;
                NewGen = newleaves;

                
            end
            
            obj.leafArray = PreviousGen;
            
            obj.TreeGrid = obj.ChebRoot.leafGrids();
            
            if ~obj.ChebRoot.is_leaf
                obj.ChebRoot.findIndex([]);
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
                
%                 WEIGHTSG = cell(1,obj.ChebRoot.dim);
%                 [WEIGHTSG{:}] = ndgrid(WEIGHTSVALS{:});
%                 
%                 WEIGHTS = ones(size(WEIGHTSG{1}));
%                 
%                 for j=1:obj.ChebRoot.dim
%                     WEIGHTS = WEIGHTS.*WEIGHTSG{j};
%                 end

                WEIGHTS = chebfun3.txm(chebfun3.txm(chebfun3.txm(ones([N,N,N]), WEIGHTSVALS{1}.', 1), ...
                    WEIGHTSVALS{2}.', 2), WEIGHTSVALS{3}.', 3);
                
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
