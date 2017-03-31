classdef PUFun < handle
    
    properties
        ChebRoot
    end
    
    methods
        
        function obj = PUFun(domain,deg_in,f)
            obj.ChebRoot = ChebPatch(domain,deg_in);
            
            %Refine on f(x)
            
            while ~obj.ChebRoot.is_refined
                
                obj.ChebRoot.sample(f);
                
                if obj.ChebRoot.is_leaf
                    obj.ChebRoot = obj.ChebRoot.splitleaf();
                else
                    obj.ChebRoot.PUsplit();
                end
                
            end
        end
        
        function ef = evalf(obj,X)
            ef = obj.ChebRoot.evalf(X);
        end
        
    end
    
end

