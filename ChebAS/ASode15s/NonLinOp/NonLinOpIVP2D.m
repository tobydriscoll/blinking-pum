classdef NonLinOpIVP2D < NonLinOp2D
    %NonLinOp Abstract class used for our PDE solvers on rectangular
    %domains.
    
    properties
        M %Mass matrix
        M_time_dep %Mass matrix is time dependent
        timederiv
        timeJac
        initial
    end
    
    methods
        
        function obj = NonLinOpIVP2D(Chebs,timederiv,timeJac,boundf,M_time_dep)
            
            @rfun
            
            obj = obj@NonLinOp2D(Chebs,timederiv,timeJac,boundf);
            obj.M_time_dep = M_time_dep;
        end
        
        function M = MassMatrix(obj,t)
            M = [];
        end
        
        % Residual used to solve for an implicit step within
        % a backwards stepping scheme (like ode15s).
        %
        % INPUT:
        %           u: solution operator will be evaluated at
        %           t: current time
        %         rhs: constant RHS of residual
        %
        % OUTPUT
        %           R: time step residual
        function R = residual(obj,u,params)
            dt = params.dt;
            if obj.M_time_dep
                MassMat = MassMatrix(obj,params.t);
                R = dt*obj.timederiv(u,params.t)-MassMat*u;
            else
                R = dt*obj.timederiv(u,params.t)-obj.M*u;
            end
        end
        
        % Calculates Jacobian for time stepping scheme.
        %
        % INPUT:
        %           u: solution operator will be evaluated at
        %           t: current time
        %        opts
        % OUTPUT
        %           R: time step residual
        function J = jac(obj,u,t,params)
            dt = params.dt;
            if obj.M_time_dep
                MassMat = MassMatrix(obj);
                J = dt*timeJac(u,t)-MassMat;
            else
                J = dt*timeJac(u,t)-obj.M;
            end
        end
    end
    
end



