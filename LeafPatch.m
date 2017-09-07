classdef LeafPatch<Patch
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties (Abstract)
        values
    end
    
    properties
        index = [];
        chebweights = [];
    end
    
    methods (Abstract)
        %This method will split the child, creating a new PUPatch. If the
        %obj does not need to split, the method returns obj.
        Child = splitleaf(obj);
        
        ef = evalf(obj,X,diff_dim,order)
    end
    
end
