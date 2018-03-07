classdef LSPatch3D < LSPatch
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
    end
    
    properties (Access = protected)
        swap_deg_in
    end
    
    methods
       % function obj = LSPatch2D(in_domain,max_lengths,domain,zone,outerbox,deg_in,cheb_deg_in,split_flag,tol,cdeg_in)
       function obj = LSPatch3D(varargin)
           
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
       end
        
        
        
        function IsGeometricallyRefined = IsLeafGeometricallyRefined(obj)
            
            lengths = [diff(obj.domain(1,:));diff(obj.domain(2,:));diff(obj.domain(3,:))];
            
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
        
        
        function max_val = sample(obj,f,grid_opt)
            
            if(nargin==2)
                grid_opt = false;
            end
            
            max_val = 0;
            
            if obj.is_geometric_refined
                
                x = chebpts(obj.degs(1)*2,obj.domain(1,:));
                y = chebpts(obj.degs(2)*2,obj.domain(2,:));
                z = chebpts(obj.degs(2)*2,obj.domain(3,:));
                
                [X,Y] = ndgrid(x,y);
                
                XP = [X(:) Y(:)];
                
                ind = obj.domain_in.Interior(XP);
                
                Mx = clenshaw(chebpts(obj.degs(1)*2),eye(obj.degs(1)));
                My = clenshaw(chebpts(obj.degs(2)*2),eye(obj.degs(2)));
                Mz = clenshaw(chebpts(obj.degs(3)*2),eye(obj.degs(3)));
                
                M = kron(Mz,kron(My,Mx));
                M = M(ind,:);
                
                if~ grid_opt
                    F = f(X(ind),Y(ind),Z(ind));
                else
                    F = f({x,y,z});
                    F = F(ind);
                end
                
                max_val = max(abs(F));
                
                warning('off','all');
                obj.coeffs = reshape(M\F,obj.degs);
                warning('on','all');
                
                E = obj.evalfGrid({x,y,z});
                E = E(ind);
                E = E(:) - F;
                
                %This is used to determin the point wise error
                obj.mid_values_err = max(abs(E(:)));  
                
            end
        end
    end
end