classdef LeafPatch<Patch
    % This is the abstract class for a leaf object. This is used with the
    % PUPatch object.
    
    properties
        index = [];
        chebweights = [];
        bump
    end
    
    methods
        % Evaluates the bump function of a patch.
        %
        %  Input:
        %      X: set of points to evaluate at
        %
        % Output:
        %     ef: array of length(X) of the patch's weight evaluated at X.
        function ef = evalfBump(obj,X)
            
            ef = ones(size(X,1),1);
            for k=1:obj.dim
                ef = ef.*obj.bump{k}(X(:,k));
            end
            
        end
        
        % Evaluates the bump function of a patch.
        %
        %  Input:
        %      X: cellarray of grids to be evaluated on.
        %
        % Output:
        %     ef: array of dim(X) of the patch's weight evaluated at X.
        function ef = evalfGridBump(obj,X)
            
            W = cell(3,1);
            
            for i=1:obj.dim
                W{i} = obj.bump{i}(X{i});
            end
            
            if obj.dim==2
                ef = W{1}*W{2}.';
            else
                ef = reshape(W{3},1,1,length(W{3})).*(W{2}'.*W{1});
            end
            
        end
    end
    methods (Abstract)
        %This method will split the child, creating a new PUPatch. If the
        %obj does not need to split, the method returns obj.
        Child = splitleaf(obj);
        
        ef = evalf(obj,X)
        
        ef = evalfGrid(obj,X)
    end
    
end
