classdef (Abstract) NonLinOp < handle
    %NonLinOp Abstract class used for our PDE solvers on rectangular
    %domains.
    
    properties
    end
    
    methods (Abstract)
        R = residual(obj,u,params);
        J = jac(obj,u,params);
    end
    
end

