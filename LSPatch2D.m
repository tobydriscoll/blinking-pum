classdef LSPatch2D<LeafPatch
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        degs %array of degrees along the dimensions
        coeffs %grid of values to be used for interpolation
        tol %tolerence used for refinement
        Out_Domain %Domain of outer square
        max_lengths %Max lengths of patch
        values
    end
    
    methods
        function obj = LSPatch2D(domain,sqr_domain,degs,max_lengths)
            obj.domain = domain;
            obj.Out_Domain = sqr_domain;
            obj.degs = degs;
            obj.cheb_length = prod(obj.degs);
            obj.is_leaf = true;
            obj.is_refined = false;
            obj.max_lengths = max_lengths;
            obj.values = [];
        end
        
        % TODO. Figure out what to do here!
        function ln=length(obj)
            ln = 0;
        end
        
        % TODO. Figure out what to do here!
        function pts = points(obj)
            pts = [];
        end
        
        % TODO. Figure out what to do here!
        function ef = evalf(obj,X,diff_dim,order)
            ef = [];
        end
        
        % TODO. Figure out what to do here!
        function ef = evalfGrid(obj,x,dim,order)
            ef = [];
        end
        
        
        function IsGeometricallyRefined = IsGeometricallyRefined(obj)
            %outer_points_s = [0.5 0.5;-0.5 0.5;0.5 -0.5;-0.5 -0.5];
            %outer_points_s = [0.75 0.75;-0.75 0.75;0.75 -0.75;-0.75 -0.75];
            %center_point = 0.5*[sum(obj.Out_Domain(1,:)) sum(obj.Out_Domain(2,:))];
            
            %outer_points(:,1) = 0.5*(diff(obj.Out_Domain(1,:))*outer_points_s(:,1)+sum(obj.Out_Domain(1,:)));
            %outer_points(:,2) = 0.5*(diff(obj.Out_Domain(2,:))*outer_points_s(:,2)+sum(obj.Out_Domain(2,:)));
            
            lengths = [diff(obj.Out_Domain(1,:));diff(obj.Out_Domain(2,:))];
            
            is_less_max = lengths<=obj.max_lengths;
            
            IsGeometricallyRefined = all(is_less_max);
            
            obj.is_geometric_refined = IsGeometricallyRefined;
        end
        
        %TODO. Factor in refinement of the function.
        function Child = splitleaf(obj)
            if obj.is_geometric_refined || IsGeometricallyRefined(obj)
                Child = obj;
            else
                %we need to split.
                lengths = [diff(obj.Out_Domain(1,:));diff(obj.Out_Domain(2,:))];
                
                is_less_max = lengths<=obj.max_lengths;
                %We split along the max length of the outer box.
                
                
                if all(is_less_max)
                    [~,split_dim] = max(lengths);
                else
                    lengths(is_less_max) = 0;
                    [~,split_dim] = max(lengths);
                end
                
                delta = 0.5*(1+PUWeights.overlap)*...
                    (obj.Out_Domain(split_dim,2)-obj.Out_Domain(split_dim,1));
                
                domain1 = obj.Out_Domain;
                domain2 = obj.Out_Domain;
                
                domain1(split_dim,:) = [obj.Out_Domain(split_dim,1) obj.Out_Domain(split_dim,1)+delta];
                domain2(split_dim,:) = [obj.Out_Domain(split_dim,2)-delta obj.Out_Domain(split_dim,2)];
                
                overlap_in = [obj.Out_Domain(split_dim,2)-delta obj.Out_Domain(split_dim,1)+delta];
                
                %We first figure out if the the subdomains sit entirely in the domain itself.
                %In this case, we would just use a standard chebyshev
                %tensor product approximation.
                x1 = chebpts(16,domain1(1,:))';
                y1 = chebpts(16,domain1(2,:))';
                [X1,Y1] = ndgrid(x1,y1);
                XP1 = [X1(:),Y1(:)];
                
                x2 = chebpts(16,domain2(1,:))';
                y2 = chebpts(16,domain2(2,:))';
                [X2,Y2] = ndgrid(x2,y2);
                XP2 = [X2(:),Y2(:)];
                
                
                
                if all(obj.domain.Interior(XP1)) 
                    %The square is in the domain. Set the child to a
                    %standard Chebpatch
                    children{1} = ChebPatch(domain1,obj.degs);
                else
                    %The square is not in the domain. Set the child to a
                    %least square patch
                    children{1} = LSPatch2D(obj.domain,domain1,obj.degs,obj.max_lengths);
                end
                
                if all(obj.domain.Interior(XP2))
                    %The square is in the domain. Set the child to a
                    %standard Chebpatch
                    children{2} = ChebPatch(domain2,obj.degs);
                else
                    %The square is not in the domain. Set the child to a
                    %least square patch
                    children{2} = LSPatch2D(obj.domain,domain2,obj.degs,obj.max_lengths);
                end
                
                x = chebpts(16,obj.Out_Domain(1,:))';
                y = chebpts(16,obj.Out_Domain(2,:))';
                [X,Y] = ndgrid(x,y);
                
                XP = [X(:),Y(:)];
                
                ind = obj.domain.Interior(XP);
                
                XP = XP(ind,:);
                
                ind11 = XP(:,split_dim)<=domain1(split_dim,2);
                
                ind22 = XP(:,split_dim)>=domain2(split_dim,1);
                
                
                outer_points_s = 0.5*[1 1;-1 1;1 -1;-1 -1];
                center_point1 = 0.5*[sum(domain1(1,:)) sum(domain1(2,:))];
                
                outer_points1(:,1) = 0.5*(diff(domain1(1,:))*outer_points_s(:,1)+sum(domain1(1,:)));
                outer_points1(:,2) = 0.5*(diff(domain1(2,:))*outer_points_s(:,2)+sum(domain1(2,:)));
                points1 = [outer_points1;center_point1];
                
                lengths1 = [diff(domain1(1,:));diff(domain1(2,:))];
                
                center_point2 = 0.5*[sum(domain2(1,:)) sum(domain2(2,:))];
                
                outer_points2(:,1) = 0.5*(diff(domain2(1,:))*outer_points_s(:,1)+sum(domain2(1,:)));
                outer_points2(:,2) = 0.5*(diff(domain2(2,:))*outer_points_s(:,2)+sum(domain2(2,:)));
                points2 = [outer_points2;center_point2];
                
                lengths2 = [diff(domain2(1,:));diff(domain2(2,:))];
                
                
                if all(ind11)
                    
                    Child = children{1};
                    
                elseif all(ind22)
                    
                    %The domain sits entirely in the second child
                    Child = children{2};
                else
                    %Return the PUPatch with the new children
                    Child = PUPatch(obj.Out_Domain,overlap_in,length(children{1})+length(children{2}),children,split_dim);
                end
            end
        end
        
        %TODO. Figure out what to do here!
        function sample(obj,f)
            
        end
        
        
        function plotdomain(obj)
            hold on;
            lengths = [diff(obj.Out_Domain(1,:));diff(obj.Out_Domain(2,:))];
            rectangle('position',[obj.Out_Domain(:,1)' lengths'],'LineWidth',2);
        end
        
        
    end
end

