classdef Disk
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radius = 1;
        center = [0 0];
    end
    
    methods
        function obj = Disk(radius,center)
            obj.radius = radius;
            obj.center = center;
        end
        
        function ind = Interior(obj,pts)
            ind = sum((pts-obj.center).^2,2)<=obj.radius.^2;
        end
        
        function plot(obj)
            
            Nxf = 66;
            Nyf = 66;
            
            x_f = chebpts(Nxf,[obj.center(1)-obj.radius obj.center(1)+obj.radius])';
            y_f = chebpts(Nyf,[obj.center(2)-obj.radius obj.center(2)+obj.radius])';
            
            [Xf,Yf] = ndgrid(x_f,y_f);
            
            XPf = [Xf(:) Yf(:)];
            
            grid_sq_ind = obj.Interior(XPf);
            
            %fine grid in domain
            XPf = XPf(grid_sq_ind,:);
            
            scatter(XPf(:,1),XPf(:,2),'red');
        end
    end
end

