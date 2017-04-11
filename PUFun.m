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
        
        function ef = evalf(obj,X,dim,order)
            ef = obj.ChebRoot.evalf(X,dim,order);
        end
                function ef = evalfGrid(obj,X,dim,order)
            ef = obj.ChebRoot.evalfGrid(X,dim,order);
        end
        
        function disp(obj)
            disp(obj.ChebRoot.toString());
        end
    end
end

