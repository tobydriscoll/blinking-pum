classdef LSPatch2D < LSPatch
    % LSPatch2D PUFun class for representing n-d functions on non-square domains.
    %
    % This class represents a single tensor product polynomial, where the
    % domain of the Chebyshev polynomial bounds the domain of the function to
    % be approximated. The coefficients are found by solving a rank deficient
    % least square problem that minimizes the l_2 norm between a given function
    % f and the Chebyshev polynomial for a set of points inside the domain of
    % f.
    
    % LSPatch2D(varargin) constructs a tensor product approximation
    % representing a function, based on the options passed into with varargin;
    % that is PUFun('perf1',perf1,'pref2',pref2,..) is called. This
    % preferences that can be set are:
    %
    % The max lengths of the patches before sampling is to occur:
    % 'MaxLengths', [d_1 d_2]
    %
    % *The non square domain: 'InnerDomain', domain object
    %
    % *The domain used for the Chebyshev polynomial: 'domain', [a,b;c,d]
    %
    % *The zone (non overlapping part from partition) used: 'zone', [a,b;c,d]
    %
    % *The domain of the root of the tree: 'outerbox', [a,b;c,d]
    %
    % *An array of boolean indicies indicating if the approximation can be
    % split in a given dimension: 'canSplit', [bool_1,bool2]
    %
    % *The tolerance used for refinement: 'tol', 1e-b
    %
    % *The degree indices from the standard degrees in each dimension for non
    % square domains : 'degreeIndex', [ind_1,ind_2].
    %
    % *The coarse degree to be used (if applicable)
    % : 'coarseDegreeIndex', [ind_1,ind_2].
    %
    % *The degree indices from the standard degrees in each dimension for
    % square domains : 'ChebDegreeIndex', [ind_1,ind_2].
    %
    % Here the degrees can be chosen from the set [3 5 9 17 33 65 129].
    % So if 'degreeIndex', [5 5 5], the max degree of any approximate will be
    % 33 in each direction.
    %
    % LSPatch2D(struct) will construct an approximation with a structure
    % struct. Here struct must contain the following fields:
    % in_domain : inner non square domain
    % outerbox : domain of the outerbox
    % zone : domain of the zone
    % domain : square domain of the polynomial
    % deg_in : indicies of degree for polynomials representing non square domains
    % cheb_deg_in : indicies of degree for polynomials representing square domains
    % cdeg_in : indicies of coarse degree
    % split_flag : boolean array indiciating to split in a dimension or not
    % max_lengths : obj.max_lengths: The max lengths of the patches before sampling is to occur
    % tol : tolerance used for refinement
    
    
    properties
        mid_values_err = inf %Store the evaluation at the Cheb points of the first kind
        Tikhonov_param = 1e-7;
        
    end
    
    properties (Access = protected)
        swap_deg_in
    end
    
    
    
    methods
        % function obj = LSPatch2D(in_domain,max_lengths,domain,zone,outerbox,deg_in,cheb_deg_in,split_flag,tol,cdeg_in)
        function obj = LSPatch2D(varargin)
            
            if length(varargin)==1
                varargin = varargin{:};
            end
            
            %Call superclass constructor
            obj = obj@LSPatch(varargin);
            
        end
        
        %Returns structure of parameters
        function p_struct = params(obj)
            
            p_struct.outerbox = obj.outerbox;
            p_struct.zone = obj.zone;
            p_struct.domain = obj.domain;
            p_struct.deg_in = obj.deg_in;
            p_struct.cheb_deg_in = obj.cheb_deg_in;
            p_struct.domain_in = obj.domain_in;
            p_struct.split_flag = obj.split_flag;
            p_struct.max_lengths = obj.max_lengths;
            p_struct.tol = obj.tol;
            p_struct.cdeg_in = obj.cdeg_in;
            p_struct.Tikhonov_param = obj.Tikhonov_param;
            
        end
        
        
        
        function Child = splitleafGeom(obj)
            
            Child = obj;
            obj.is_geometric_refined = true;
            
               x = chebpts(32,obj.zone(1,:))';
               y = chebpts(32,obj.zone(2,:))';
               
               [X,Y] = ndgrid(x,y);
                        
               XP = [X(:) Y(:)];
                        
               ind = obj.domain_in.Interior(XP);
               
               XP = XP(ind,:);
               
               new_zone = zeros(2,2);
               
               new_zone(:,1) = min(XP);
               
               new_zone(:,2) = max(XP);
               
               %pudge out zone a bit
               deltax = 0.25*obj.overlap*diff(obj.zone(1,:));
                %The width of the overlap
               deltay = 0.25*obj.overlap*diff(obj.zone(2,:));       
               
               new_zone(1,1) = max(new_zone(1,1)-deltax,obj.zone(1,1));
               new_zone(1,2) = min(new_zone(1,2)+deltax,obj.zone(1,2));
               
               new_zone(2,1) = max(new_zone(2,1)-deltay,obj.zone(2,1));
               new_zone(2,2) = min(new_zone(2,2)+deltay,obj.zone(2,2));
               
               %The width of the overlap
               deltax = 0.25*obj.overlap*diff(obj.zone(1,:));
                %The width of the overlap
               deltay = 0.25*obj.overlap*diff(obj.zone(2,:));

               obj.zone = new_zone;
               
               new_zone(1,1) = max(new_zone(1,1)-deltax,obj.outerbox(1,1));
               new_zone(1,2) = min(new_zone(1,2)+deltax,obj.outerbox(1,2));
               
               new_zone(2,1) = max(new_zone(2,1)-deltay,obj.outerbox(2,1));
               new_zone(2,2) = min(new_zone(2,2)+deltay,obj.outerbox(2,2));
               
               obj.domain = new_zone;

               
%                 lengths = [diff(obj.domain(1,:));diff(obj.domain(2,:))];
%                 
%                 is_less_max = lengths<=obj.max_lengths;
%                 
%                 is_less_max = all(is_less_max);
%                 
%                 in_sub = false(4,1);
%                 
%                 if is_less_max
%                     
%                     splitx = mean(obj.domain(1,:));
%                     splity = mean(obj.domain(2,:));
%                     
%                     subdomain{1} = [obj.domain(1,1) splitx;obj.domain(2,1) splity];
%                     subdomain{2} = [obj.domain(1,1) splitx;splity obj.domain(2,2)];
%                     subdomain{3} = [splitx obj.domain(1,2);obj.domain(2,1) splity];
%                     subdomain{4} = [splitx obj.domain(1,2);splity obj.domain(2,2)];
%                     
%                     for i=1:4
%                         c_subdom = subdomain{i};
%                         
%                         x = chebpts(16,c_subdom(1,:))';
%                         y = chebpts(16,c_subdom(2,:))';
%                         
%                         [X,Y] = ndgrid(x,y);
%                         
%                         XP = [X(:) Y(:)];
%                         
%                         in_sub(i) = sum(obj.domain_in.Interior(XP))>40;
%                     end
%                 end
%                 
%                 obj.is_geometric_refined = is_less_max & sum(in_sub)>=2;
%                     
%                 Child = obj;
%                 if ~obj.is_geometric_refined
%                     for k=1:obj.dim
%                         if Child.is_leaf
%                             Child = split(Child,k);
%                         else
%                             Child.split(k);
%                         end
%                     end
%                 end 
        end
        
        function Child = splitleaf(obj,Max,set_vals)
            
            obj.GlobalMax = Max;
            
            if obj.mid_values_err>obj.tol
                
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
                %Chop(obj);
                Child = obj;
                Child.is_refined = true;
            end
            
        end
        
        
        
        function max_val = sample(obj,f,grid_opt)
            
            
            
            if(nargin==2)
                grid_opt = false;
            end
            
            max_val = 0;
            
            if obj.is_geometric_refined
                
                x = chebpts(obj.degs(1)*2,obj.domain(1,:));
                y = chebpts(obj.degs(2)*2,obj.domain(2,:));
                
                [X,Y] = ndgrid(x,y);
                
                XP = [X(:) Y(:)];
                
                ind = obj.domain_in.Interior(XP);
                
%                 D_x = ChebDiff(obj.degs(1));
%                 D_y = ChebDiff(obj.degs(2));
%                 D_xx = ChebDiff(obj.degs(1))^2;
%                 D_yy = ChebDiff(obj.degs(2))^2;
                
               % Lap = kron(D_yy,eye(obj.degs(1)))+kron(eye(obj.degs(2)),D_xx)+2*kron(D_y,D_x);
                
                Mx = clenshaw(chebpts(obj.degs(1)*2),eye(obj.degs(1)));
                My = clenshaw(chebpts(obj.degs(2)*2),eye(obj.degs(2)));
                
                M = kron(My,Mx);
                M = M(ind,:);
                
                if~ grid_opt
                    F = f(X(ind),Y(ind));
                else
                    F = f({x,y});
                    F = F(ind);
                end
                
                max_val = max(abs(F));
                
                obj.LocalMax = max_val;
                
                %M2 = [obj.Tikhonov_param*Lap;M];
                %F2 = [zeros(prod(obj.degs),1);F];
                %warning('off','all');
                %obj.coeffs = reshape(M2\F2,obj.degs);
                %warning('on','all');
                
                 warning('off','all');
                 obj.coeffs = reshape(M\F,obj.degs);
                 warning('on','all');
                
                E = obj.evalfGrid({x,y});
                E = E(ind);
                E = E(:) - F;
                
                %This is used to determin the point wise error
                obj.mid_values_err = max(abs(E(:)));
                
            end
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
            
            struct0 = obj.params;
            struct0.domain = domain0; struct0.zone = zone0;
            
            struct1 = obj.params;
            struct1.domain = domain1; struct1.zone = zone1;
            
            if all(obj.domain_in.Interior(XP1))
                %The square is in the domain. Set the child to a
                %standard Chebpatch
                children{1} = ChebPatch(struct0);
            else
                %The square is not in the domain. Set the child to a
                %least square patch
                children{1} = LSPatch2D(struct0);
                children{1}.is_geometric_refined = false;
            end
            
            if all(obj.domain_in.Interior(XP2))
                %The square is in the domain. Set the child to a
                %standard Chebpatch
                children{2} = ChebPatch(struct1);
            else
                %The square is not in the domain. Set the child to a
                %least square patch
                children{2} = LSPatch2D(struct1);
                children{2}.is_geometric_refined = false;
            end
            
            x = chebpts(16,obj.domain(1,:))';
            y = chebpts(16,obj.domain(2,:))';
            
            [X,Y] = ndgrid(x,y);
            
            XP = [X(:),Y(:)];
            
            ind = obj.domain_in.Interior(XP);
            
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
                Child = PUPatch(obj.domain,obj.zone,children,split_dim);
                
                Child.is_geometric_refined = false;
            end
            

            if set_vals
                for k=1:2
                    Child.children{k}.sample(obj.evalfGrid(Child.children{k}.leafGrids()));
                end
            end
        end
        
        function Chop(obj)
            
            loc_tol = obj.tol^(7/8);
            
            for k=1:obj.dim
                
                if obj.split_flag(k)
                    
                    colChebtech = chebfun3t.unfold(obj.coeffs, k);
                    colChebtech = sum(abs(colChebtech),2);
                    fCol = chebtech2({[],colChebtech});
                    hscale = diff(obj.domain(k,:));
                    
                    tol = loc_tol*max(obj.GlobalMax/obj.LocalMax,hscale);
                    cutoff = length(simplify(fCol, tol))+1;
                    
                    if cutoff<obj.degs(k)
                        j = find(cutoff<=obj.standard_degs);
                        
                        if j<obj.deg_in(k)
                            obj.deg_in(k) = j;
                            obj.deg(k) = obj.standard_degs(j);
                            
                            if k==1
                                obj.coeffs = obj.coeffs(1:obj.deg(k),:);
                            else
                                obj.coeffs = obj.coeffs(:,1:obj.deg(k));
                            end
                        end
                    end
                end
                
            end
            
        end
        
    end
end