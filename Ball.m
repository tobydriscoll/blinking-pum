classdef Ball
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radius = 1;
        center = [0 0 0];
    end
    
    methods
        function obj = Ball(radius,center)
            obj.radius = radius;
            obj.center = center;
        end
        
        function ind = Interior(obj,pts)
            ind = sum((pts-obj.center).^2,2)<=obj.radius.^2;
        end
        
        function plot(obj)
            
        end
    end
end

