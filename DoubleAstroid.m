classdef DoubleAstroid
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = DoubleAstroid()
        end
        
        function ind = Interior(obj,pts)
            
            [THo,R] = cart2pol(pts(:,1),pts(:,2));
            
            TH = mod(THo-pi+pi/8,2*pi)-pi;
            
            TH = mod(TH,pi/2);
            
            ind1 = R<=(cos(TH).^(2/3)+sin(TH).^(2/3)).^(-3/2);
           
            TH = mod(THo-pi-pi/8,2*pi)-pi;
            
            TH = mod(TH,pi/2);
            
            ind2 = R<=(cos(TH).^(2/3)+sin(TH).^(2/3)).^(-3/2);
            
            ind = ind1 | ind2;
            
        end
        
        function plot(obj)
            
            Nxf = 240;
            Nyf = 240;
            
            x_f = linspace(-1,1,Nxf)';
            y_f = linspace(-1,1,Nyf)';
            
            [Xf,Yf] = ndgrid(x_f,y_f);
            
            XPf = [Xf(:) Yf(:)];
            
            grid_sq_ind = obj.Interior(XPf);
            
            %fine grid in domain
            XPf = XPf(grid_sq_ind,:);
            
            scatter(XPf(:,1),XPf(:,2),'red');
        end
    end
end

