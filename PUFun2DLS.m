classdef PUFun2DLS < handle
    
    properties
        ChebRoot
        TreeGrid
    end
    
    methods
        
        function obj = PUFun2DLS(domain,sqr_domain,degs,max_lengths)
            obj.ChebRoot = LSPatch2D(domain,sqr_domain,degs,max_lengths);
            
            while ~obj.ChebRoot.is_geometric_refined
                
                if obj.ChebRoot.is_leaf
                    obj.ChebRoot = obj.ChebRoot.splitleaf();
                else
                    obj.ChebRoot.PUsplit();
                end
                
            end
        end
        
    end
end