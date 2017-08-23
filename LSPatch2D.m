classdef LSPatch2D<LeafPatch
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        degs %array of degrees along the dimensions
        coeffs %grid of values to be used for interpolation
        Out_Domain %Domain of outer square
        max_lengths %Max lengths of patch
        values
        pinvM %Store the pseudoInverse
        mid_values_err = inf %Store the evaluation at the Cheb points of the first kind
        cheblength = 32;
    end
    
    properties (Constant)
        invf = @(x,dom) 2/diff(dom)*x-sum(dom)/diff(dom); %takes points from a domain to [-1 1]
        forf = @(x,dom) 0.5*diff(dom)*x+0.5*sum(dom); %takes points from [-1 1] to a domain
    end
    
    methods
        function obj = LSPatch2D(domain,sqr_domain,degs,max_lengths)
            obj.domain = domain;
            obj.Out_Domain = sqr_domain;
            obj.degs = degs;
            obj.cheb_length = prod(obj.degs);
            obj.is_leaf = true;
            obj.max_lengths = max_lengths;
            obj.tol = 1e-12;
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
        
        function ef = evalfGrid(obj,X,diff_dim,order)
            
            G = obj.coeffs;
            
            for k=1:2
                %Shift the points to the domain [-1 1]x[-1 1]
                X{k} = obj.invf(X{k},obj.Out_Domain(:,k));
                
                %Evaluate the points at the Chebyshev polynomials
                F = clenshaw(X{k},eye(obj.cheblength));
                
                %Multiply the coefficients with F
                G = chebfun3t.txm(G, F, k);
            end
            
            ef  = G;
            
        end
        
        
        function IsGeometricallyRefined = IsGeometricallyRefined(obj)
            outer_points_s = [1 1;-1 1;1 -1;-1 -1];
            
            outer_points(:,1) = obj.forf(outer_points_s(:,1),obj.Out_Domain(1,:));
            outer_points(:,2) = obj.forf(outer_points_s(:,2),obj.Out_Domain(2,:));
            
            lengths = [diff(obj.Out_Domain(1,:));diff(obj.Out_Domain(2,:))];
            
            is_less_max = lengths<=obj.max_lengths;
            
            %IsGeometricallyRefined = all(is_less_max) && any(obj.domain.Interior(outer_points));
            IsGeometricallyRefined = all(is_less_max);
            obj.is_geometric_refined = IsGeometricallyRefined;
        end
        
        function Child = splitleaf(obj)
            
            if ~obj.is_geometric_refined
                obj.IsGeometricallyRefined();
            end
            
            
            obj.is_refined = obj.mid_values_err<= obj.tol;
            %obj.is_refined = true;
            
            if obj.is_geometric_refined
                a=1;
            end
            
            
            if obj.is_geometric_refined && obj.is_refined
                
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
                
                
                if all(ind11)
                     %The domain sits entirely in the first child
                    Child = children{1};
                    
                elseif all(ind22)
                    
                    %The domain sits entirely in the second child
                    Child = children{2};
                else
                    %Return the PUPatch with the new children
                    Child = PUPatch(obj.Out_Domain,overlap_in,length(children{1})+length(children{2}),children,split_dim,obj.index);
                end
            end
        end
        
        function sample(obj,f)
            if obj.is_geometric_refined
                
                x = chebpts(obj.cheblength*2,obj.Out_Domain(1,:));
                y = chebpts(obj.cheblength*2,obj.Out_Domain(2,:));
                
                x1 = chebpts(obj.cheblength*2,obj.Out_Domain(1,:),1);
                y1 = chebpts(obj.cheblength*2,obj.Out_Domain(2,:),1);
                
                [X,Y] = ndgrid(x,y);
                
                [X1,Y1] = ndgrid(x1,y1);
                
                XP = [X(:) Y(:)];
                
                XP1 = [X1(:) Y1(:)];
                
                ind = obj.domain.Interior(XP);
                ind1 = obj.domain.Interior(XP1);
                
                XP = XP(ind,:);
                
                if ~obj.is_refined   
                    Mx = clenshaw(chebpts(obj.cheblength*2),eye(obj.cheblength));
                    M = kron(Mx,Mx);
                    obj.pinvM = pinv(M(ind,:));
                end
                
                obj.coeffs = reshape(obj.pinvM*f(XP),[obj.cheblength obj.cheblength]);
                
                E = obj.evalfGrid({x1,y1},1,0);
                E = E(:) - f(XP1);
                E = E(ind1);
                
                %This is used to determin
                obj.mid_values_err = max(abs(E(:)));
            end
        end
        
        
        function plotdomain(obj)
            hold on;
            lengths = [diff(obj.Out_Domain(1,:));diff(obj.Out_Domain(2,:))];
            rectangle('position',[obj.Out_Domain(:,1)' lengths'],'LineWidth',2);
            hold off;
        end
        
        
    end
end

