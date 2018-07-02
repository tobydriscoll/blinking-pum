classdef Diamond
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = Diamond()
        end
        
        function ind = Interior(obj,pts)
            
            ind = abs(pts(:,1))+abs(pts(:,2))<1;
            
        end
        
        function bound = Boundary(obj,N)
            
            num_b = floor(N/4);
            x = linspace(-sqrt(2)/2,sqrt(2)/2,num_b)';
            
            xup = flipud(x);
            
            P = [x -sqrt(2)/2*ones(num_b,1);sqrt(2)/2*ones(num_b-1,1) x(2:end);xup(2:end) sqrt(2)/2*ones(num_b-1,1);-sqrt(2)/2*ones(num_b-2,1) xup(2:end-1)];
            
            R = [cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)];
            
            bound = (R*P')';
            
        end
        
        function plot(obj)
            
            Nxf = 200;
            Nyf = 200;
            
            x_f = linspace(-1,1,Nxf)';
            y_f = linspace(-1,1,Nyf)';
            
            [Xf,Yf] = ndgrid(x_f,y_f);
            
            ind = obj.Interior([Xf(:) Yf(:)]);
            
            P = [Xf(ind) Yf(ind)];
            
            TRI = delaunay(P);
            
            Z = zeros(length(P),1);
            
            trisurf(TRI,Xf(ind),Yf(ind),Z);
            view(0,90);
            shading interp;
        end
    end
end

