classdef LSPatch2D<LeafPatch
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        degs %array of degrees along the dimensions
        coeffs %grid of values to be used for interpolation
        in_domain
        max_lengths %Max lengths of patch
        values
        pinvM %Store the pseudoInverse
        mid_values_err = inf %Store the evaluation at the Cheb points of the first kind
        deg_in
        cheblength = 17;
    end
    
    properties (Access = protected)
        cdeg_in %index for the course degrees
        swap_deg_in
    end
    
    properties (Constant)
        standard_variables = load('cheb_points_matrices.mat');
        standard_degs = [3 5 9 17 33 65 129];
        invf = @(x,dom) 2/diff(dom)*x-sum(dom)/diff(dom); %takes points from a domain to [-1 1]
        forf = @(x,dom) 0.5*diff(dom)*x+0.5*sum(dom); %takes points from [-1 1] to a domain
    end
    
    methods
        function obj = LSPatch2D(in_domain,max_lengths,domain,zone,outerbox,deg_in,split_flag,tol,cdeg_in)
            
            %Call superclass constructor
            obj = obj@LeafPatch(domain,zone,outerbox);
            obj.in_domain = in_domain;
            obj.max_lengths = max_lengths;
            
            if nargin < 6
                obj.deg_in = zeros(1,obj.dim);
                obj.cdeg_in = zeros(1,obj.dim);
                obj.deg_in(:) = 7;
                obj.cdeg_in(:) = 3;
                obj.split_flag = true(obj.dim,1);
                obj.tol = 1e-12;
            elseif nargin < 7
                obj.deg_in = deg_in;
                obj.cdeg_in = zeros(1,obj.dim);
                obj.cdeg_in(:) = 3;
                obj.split_flag = true(obj.dim,1);
                obj.tol = 1e-12;
            elseif nargin < 8
                obj.deg_in = deg_in;
                obj.cdeg_in = zeros(1,obj.dim);
                obj.cdeg_in(:) = 3;
                obj.split_flag = split_flag;
                obj.tol = 1e-12;
            elseif nargin < 9
                obj.deg_in = deg_in;
                obj.cdeg_in = zeros(1,obj.dim);
                obj.cdeg_in(:) = 3;
                obj.split_flag = split_flag;
                obj.tol = tol;
            elseif nargin < 10
                obj.deg_in = deg_in;
                obj.cdeg_in = cdeg_in;
                obj.split_flag = split_flag;
                obj.tol = tol;
            end
            
            obj.degs = obj.standard_degs(obj.deg_in);
            
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
                X{k} = obj.invf(X{k},obj.domain(k,:));
                
                %Evaluate the points at the Chebyshev polynomials
                F = clenshaw(X{k},eye(obj.degs(k)));
                
                %Multiply the coefficients with F
                G = chebfun3t.txm(G, F, k);
            end
            
            ef  = G;
            
        end
        
        
        function IsGeometricallyRefined = IsLeafGeometricallyRefined(obj)
            
            lengths = [diff(obj.domain(1,:));diff(obj.domain(2,:))];
            
            is_less_max = lengths<=obj.max_lengths;
            
            %IsGeometricallyRefined = all(is_less_max) && any(obj.domain.Interior(outer_points));
            IsGeometricallyRefined = all(is_less_max);
            obj.is_geometric_refined = IsGeometricallyRefined;
        end
        
        function Child = splitleaf(obj,Max,set_vals)
            
            if ~obj.is_geometric_refined || obj.mid_values_err>obj.tol
                
                Child = obj;
                %Go through and split in each unresolved direction
                for k=1:obj.dim
                    if Child.is_leaf
                        Child = split(Child,k);
                    else
                        Child.split(k);
                    end
                end
                
            else
                Child = obj;
                Child.is_refined = true;
            end
            
%             vscale = Max;
%             
%             loc_tol = obj.tol^(7/8);
%             
%             local_max = Max;
%             
%             
%             for k=1:obj.dim
%                 
%                 if obj.split_flag(k)
%                     
%                     colChebtech = chebfun3t.unfold(obj.coeffs, k);
%                     colChebtech = sum(abs(colChebtech),2);
%                     fCol = chebtech2({[],colChebtech});
%                     hscale = diff(obj.domain(k,:));
%                     
%                     tol = loc_tol*max(vscale./local_max,hscale);
%                     cutoff = length(simplify(fCol, tol))+1;
%                     
%                     if cutoff<obj.degs(k)
%                         obj.split_flag(k) = false;
%                         j = find(cutoff<=obj.standard_degs,1);
%                         obj.deg_in(k) = j;
%                         obj.degs(k) = obj.standard_degs(j);
%                         
%                         if k==1
%                             obj.coeffs = obj.coeffs(1:obj.degs(k),:);
%                         else
%                             obj.coeffs = obj.coeffs(:,1:obj.degs(k));
%                         end
%                     end
%                 end
%                 
%             end
%             
%             
%             
%             
%             
%             if ~any(obj.split_flag)
%                 %The leaf is refined, so return it.
%                 obj.is_refined = true;
%                 obj.cheb_length = prod(obj.degs);
%                 Child = obj;
%             else
%                 
%                 %                 ind = 1:obj.dim;
%                 %                 ind = ind(obj.split_flag);
%                 %
%                 %                 [~,split_dim] = max(diff(obj.domain(obj.split_flag,:).',1));
%                 %                 split_dim = ind(split_dim);
%                 %
%                 %                 Child = split(obj,split_dim,set_vals);
%                 
%                 
%                 Child = obj;
%                 %Go through and split in each unresolved direction
%                 for k=1:obj.dim
%                     if obj.split_flag(k)
%                         if Child.is_leaf
%                             Child = split(Child,k);
%                         else
%                             Child.split(k);
%                         end
%                     end
%                 end
%                 
%                 
%                 
%             end
        end
        
                % The method determines will split a child into along
        % a dimension.
        %
        %     Input:
        %   overlap: overlap intended to be used for the splitting
        %
        %    Output:
        %     Child: the PUPatch with the two new children.
        function Child = split(obj,split_dim,set_vals)
            
            if nargin == 2
                set_vals = false;
            end
            
            children = cell(1,2);
            
            %The width of the overlap
            delta = 0.5*obj.overlap*diff(obj.zone(split_dim,:));
            
            zone0 = obj.zone;
            zone1 = obj.zone;
            
            domain0 = obj.domain;
            domain1 = obj.domain;
            
            m = sum(obj.zone(split_dim,:))/2;
            
            zone0(split_dim,:) = [obj.zone(split_dim,1) m];
            zone1(split_dim,:) = [m obj.zone(split_dim,2)];
            
            domain0(split_dim,:) = [max(obj.outerbox(split_dim,1),obj.zone(split_dim,1)-delta) m+delta];
            domain1(split_dim,:) = [m-delta,min(obj.outerbox(split_dim,2),obj.zone(split_dim,2)+delta)];
            
            %We first figure out if the the subdomains sit entirely in the domain itself.
            %In this case, we would just use a standard chebyshev
            %tensor product approximation.
            x1 = chebpts(16,domain0(1,:))';
            y1 = chebpts(16,domain0(2,:))';
            [X1,Y1] = ndgrid(x1,y1);
            XP1 = [X1(:),Y1(:)];
            
            x2 = chebpts(16,domain1(1,:))';
            y2 = chebpts(16,domain1(2,:))';
            [X2,Y2] = ndgrid(x2,y2);
            XP2 = [X2(:),Y2(:)];
            
            if all(obj.in_domain.Interior(XP1))
                %The square is in the domain. Set the child to a
                %standard Chebpatch
                children{1} = ChebPatch(domain0,zone0,obj.outerbox,obj.deg_in,true(obj.dim,1),obj.tol,obj.cdeg_in);
            else
                %The square is not in the domain. Set the child to a
                %least square patch
                children{1} = LSPatch2D(obj.in_domain,obj.max_lengths,domain0,zone0,obj.outerbox,obj.deg_in,obj.split_flag,obj.tol,obj.cdeg_in);
            end
            
            if all(obj.in_domain.Interior(XP2))
                %The square is in the domain. Set the child to a
                %standard Chebpatch
                children{2} = ChebPatch(domain1,zone1,obj.outerbox,obj.deg_in,true(obj.dim,1),obj.tol,obj.cdeg_in);
            else
                %The square is not in the domain. Set the child to a
                %least square patch
                children{2} = LSPatch2D(obj.in_domain,obj.max_lengths,domain1,zone1,obj.outerbox,obj.deg_in,obj.split_flag,obj.tol,obj.cdeg_in);
            end
            
            x = chebpts(16,obj.domain(1,:))';
            y = chebpts(16,obj.domain(2,:))';
            
            [X,Y] = ndgrid(x,y);
            
            XP = [X(:),Y(:)];
            
            ind = obj.in_domain.Interior(XP);
            
            XP = XP(ind,:);
            
            ind11 = XP(:,split_dim)<=domain0(split_dim,2);
            ind22 = XP(:,split_dim)>=domain1(split_dim,1);
            
            
            if all(ind11)
                %The domain sits entirely in the first child
                Child = children{1};
                
            elseif all(ind22)
                
                %The domain sits entirely in the second child
                Child = children{2};
            else
                %Return the PUPatch with the new children
                Child = PUPatch(obj.domain,obj.zone,length(children{1})+length(children{2}),children,split_dim,obj.index);
            end
            
            if set_vals
                for k=1:2
                    Child.children{k}.sample(obj.evalfGrid(Child.children{k}.leafGrids()));
                end
            end
        end
        
        
        function max_val = sample(obj,f,grid_opt)
            
            if(nargin==2)
                grid_opt = false;
            end
            
            max_val = inf;
            
            if obj.is_geometric_refined
                
                x = chebpts(obj.degs(1)*2,obj.domain(1,:));
                y = chebpts(obj.degs(2)*2,obj.domain(2,:));
                
                x1 = chebpts(obj.degs(1)*2,obj.domain(1,:),1);
                y1 = chebpts(obj.degs(2)*2,obj.domain(2,:),1);
                
                [X,Y] = ndgrid(x,y);
                
                [X1,Y1] = ndgrid(x1,y1);
                
                XP = [X(:) Y(:)];
                
                XP1 = [X1(:) Y1(:)];
                
                ind = obj.in_domain.Interior(XP);
                ind1 = obj.in_domain.Interior(XP1);
                
                XP = XP(ind,:);
                
                %                 if ~obj.is_refined
                %                     Mx = clenshaw(chebpts(obj.cheblength*2),eye(obj.cheblength));
                %                     M = kron(Mx,Mx);
                %                     obj.pinvM = pinv(M(ind,:));
                %                 end
                
                Mx = clenshaw(chebpts(obj.degs(1)*2),eye(obj.degs(1)));
                My = clenshaw(chebpts(obj.degs(2)*2),eye(obj.degs(2)));
                
                M = kron(My,Mx);
                M = M(ind,:);
                
                %                 P = M*M'-eye(length(XP));
                %                 y = pinv(P*M)*P*f(XP);
                %                 z = M'*(f(XP)-M*y);
                %                 obj.coeffs = reshape(y+z,[obj.cheblength obj.cheblength]);
                
                %                obj.coeffs = reshape(pinv(M)*f(XP),[obj.cheblength obj.cheblength]);
                warning('off','all');
                obj.coeffs = reshape(M\f(XP),[obj.degs(1) obj.degs(2)]);
                warning('on','all');
                
                E = obj.evalfGrid({x1,y1},1,0);
                E = E(:) - f(XP1);
                E = E(ind1);
                
                %This is used to determin the point wise error
                obj.mid_values_err = max(abs(E(:)));
                
                
            end
        end
    end
end