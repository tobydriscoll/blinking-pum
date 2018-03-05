classdef LSPatch2D<LeafPatch
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
        degs %array of degrees along the dimensions
        cdegs
        coeffs %grid of values to be used for interpolation
        in_domain
        max_lengths %Max lengths of patch
        values
        mid_values_err = inf %Store the evaluation at the Cheb points of the first kind
        deg_in
        cheb_deg_in
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
       % function obj = LSPatch2D(in_domain,max_lengths,domain,zone,outerbox,deg_in,cheb_deg_in,split_flag,tol,cdeg_in)
       function obj = LSPatch2D(varargin)
           
           if length(varargin)==1
               varargin = varargin{:};
           end
           
           %Call superclass constructor
           obj = obj@LeafPatch(varargin);
           
            
            obj.tol = 1e-6;
            obj.deg_in = zeros(obj.dim,1);
            obj.deg_in(:) = 4;
            
            obj.cheb_deg_in = zeros(obj.dim,1);
            obj.cheb_deg_in(:) = 7;
            
            obj.cdeg_in = zeros(obj.dim,1);
            obj.cdeg_in(:) = 3;
            obj.split_flag = true(obj.dim,1);
            obj.max_lengths = inf(obj.dim,1);
            
           if isstruct(varargin)

               obj.deg_in = varargin.deg_in;
               obj.cheb_deg_in = varargin.cheb_deg_in;
               obj.in_domain = varargin.in_domain;
               obj.split_flag = varargin.split_flag;
               obj.max_lengths = varargin.max_lengths;
               obj.tol = varargin.tol;
               obj.cdeg_in = varargin.cdeg_in;
               
           else
               
               args = varargin;
               
               while ( ~isempty(args) )
                   if strcmpi(args{1}, 'degreeIndex')
                       if numel(args{2})==1
                           obj.degs_in = zeros(1,obj.dim);
                           obj.deg_in(:) = args{2};
                       else
                           obj.deg_in = args{2};
                       end
                   elseif strcmpi(args{1}, 'canSplit')
                       obj.split_flag = args{2};
                   elseif strcmpi(args{1}, 'tol')
                       obj.tol = args{2};
                   elseif strcmpi(args{1}, 'coarseDegreeIndex')
                       if numel(args{2})==1
                           obj.cdegs_in = zeros(1,obj.dim);
                       else
                           obj.cdegs_in = args{2};
                           obj.cdeg_in(:) = args{2};
                       end
                   elseif strcmpi(args{1}, 'MaxLengths')
                       if numel(args{2})==1
                           obj.max_lengths = zeros(1,obj.dim);
                           obj.max_lengths(:) = args{2};
                       else
                           obj.max_lengths = args{2};
                       end
                   elseif strcmpi(args{1}, 'ChebDegreeIndex')
                       if numel(args{2})==1
                           obj.cheb_deg_in = zeros(1,obj.dim);
                           obj.cheb_deg_in(:) = args{2};
                       else
                           obj.cheb_deg_in = args{2};
                       end
                   elseif strcmpi(args{1}, 'InnerDomain')
                       obj.in_domain = args{2};
                   elseif strcmpi(args{1}, 'domain')
                       obj.domain = args{2};
                   end
                   args(1:2) = [];
               end
               
           end
           
           
           obj.degs = obj.standard_degs(obj.deg_in);
           obj.cdegs = obj.standard_degs(obj.cdeg_in);
           
           obj.cheb_length = prod(obj.degs);
           
       end
        
       %Returns structure of parameters
       function p_struct = params(obj) 
           
           p_struct.outerbox = obj.outerbox;
           p_struct.zone = obj.zone;
           p_struct.domain = obj.domain;
           p_struct.deg_in = obj.deg_in;
           p_struct.cheb_deg_in = obj.cheb_deg_in;
           p_struct.in_domain = obj.in_domain;
           p_struct.split_flag = obj.split_flag;
           p_struct.max_lengths = obj.max_lengths;
           p_struct.tol = obj.tol;
           p_struct.cdeg_in = obj.cdeg_in;
       end
        
        % TODO. Figure out what to do here!
        function ln=length(obj)
            ln = prod(obj.degs);
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
            
            if all(obj.in_domain.Interior(XP1))
                %The square is in the domain. Set the child to a
                %standard Chebpatch
                children{1} = ChebPatch(struct0);
            else
                %The square is not in the domain. Set the child to a
                %least square patch
                children{1} = LSPatch2D(struct0);
            end
            
            if all(obj.in_domain.Interior(XP2))
                %The square is in the domain. Set the child to a
                %standard Chebpatch
                children{2} = ChebPatch(struct1);
            else
                %The square is not in the domain. Set the child to a
                %least square patch
                children{2} = LSPatch2D(struct1);
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
                Child = PUPatch(obj.domain,obj.zone,children,split_dim);
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
            
            max_val = 0;
            
            if obj.is_geometric_refined
                
                x = chebpts(obj.degs(1)*2,obj.domain(1,:));
                y = chebpts(obj.degs(2)*2,obj.domain(2,:));

                [X,Y] = ndgrid(x,y);
                
                XP = [X(:) Y(:)];
                
                ind = obj.in_domain.Interior(XP);
                
                XP = XP(ind,:);
                
                Mx = clenshaw(chebpts(obj.degs(1)*2),eye(obj.degs(1)));
                My = clenshaw(chebpts(obj.degs(2)*2),eye(obj.degs(2)));
                
                M = kron(My,Mx);
                M = M(ind,:);
                
                warning('off','all');

                XPC = num2cell(XP,1);
                F = f(XPC{:});
                
                max_val = max(abs(F));
                
                obj.coeffs = reshape(M\F,[obj.degs(1) obj.degs(2)]);
                
                %obj.coeffs = reshape(ResitrictedLS(M,F,eye(prod(obj.degs))),[obj.degs(1) obj.degs(2)]);

                warning('on','all');
                
                E = obj.evalfGrid({x,y},1,0);
                E = E(ind);
                E = E(:) - F;
                
                %This is used to determin the point wise error
                obj.mid_values_err = max(abs(E(:)));
                
                
            end
        end
    end
end