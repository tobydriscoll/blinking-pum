classdef Astroid
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        angle
    end
    
    methods
        function obj = Astroid(angle)
            obj.angle = angle;
        end
        
        function ind = Interior(obj,pts)
            
            %ang = obj.angle;
            %ROT = [cos(ang) -sin(ang);sin(ang) cos(ang)];
            
            %pts = (ROT*pts.').';
            
            [TH,R] = cart2pol(pts(:,1),pts(:,2));
            
            TH = mod(TH-pi+obj.angle,2*pi)-pi;
            
            TH = mod(TH,pi/2);
            
            ind = R<=(cos(TH).^(2/3)+sin(TH).^(2/3)).^(-3/2);
       
%             TH1 =  TH>=-pi & TH<=-pi/2;
%             TH2 =  TH>=-pi/2 & TH<=0;
%             TH3 =  TH>=0 & TH<=pi/2;
%             TH4 =  TH>=pi/2;
            
            
           
%             R1 = (cos(TH(TH1)+pi).^(2/3)+sin(TH(TH1)+pi).^(2/3)).^(-3/2);
%             R2 = (cos(TH(TH2)+pi/2).^(2/3)+sin(TH(TH2)+pi/2).^(2/3)).^(-3/2);
%             R3 = (cos(TH(TH3)).^(2/3)+sin(TH(TH3)).^(2/3)).^(-3/2);
%             R4 = (cos(TH(TH4)-pi/2).^(2/3)+sin(TH(TH4)-pi/2).^(2/3)).^(-3/2);
%             
%             ind = true(size(pts,1),1);
%             
%             ind(TH1) = R(TH1)<=R1;
%             ind(TH2) = ind(TH2) & R(TH2)<=R2;
%             ind(TH3) = ind(TH3) & R(TH3)<=R3;
%             ind(TH4) = ind(TH4) & R(TH4)<=R4;
            
            
        end
        
        function plot(obj)
            
            Nxf = 200;
            Nyf = 200;
            
            x_f = chebpts(Nxf,[-1 1])';
            y_f = chebpts(Nyf,[-1 1])';
            
            [Xf,Yf] = ndgrid(x_f,y_f);
            
            XPf = [Xf(:) Yf(:)];
            
            grid_sq_ind = obj.Interior(XPf);
            
            %fine grid in domain
            XPf = XPf(grid_sq_ind,:);
            
            scatter(XPf(:,1),XPf(:,2),'red');
        end
    end
end

