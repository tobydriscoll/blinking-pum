classdef LSPatch < LeafPatch
% LSPatch2D PUFun class for representing n-d functions on non-square domains.
%
% This class represents a single tensor product polynomial, where the
% domain of the Chebyshev polynomial bounds the domain of the function to
% be approximated. The coefficients are found by solving a rank deficient 
% least square problem that minimizes the l_2 norm between a given function
% f and the Chebyshev polynomial for a set of points inside the domain of
% f.
 
% LSPatch2D(var_struct) constructs a tensor product approximation
% representing a function, based on the options passed with the structure
% var_struct. The options are:
% 
% *var_struct.max_lengths, the max lengths of the patches before sampling 
% is to occur:.
%
% *var_struct.domain_in, the non square domain
%
% *var_struct.domain, the domain used for the Chebyshev polynomial.
%
% *var_struct.zone ,the zone (non overlapping part from partition) used.
%
% *var_struct.outerbox, the domain of the root of the tree: 'outerbox'.
%
% *var_struct.split_flag, an array of boolean indicies indicating if the 
% approximation can be split in a given dimension.
%
% *var_struct.tol, the tolerance used for refinement: 'tol', 1e-b
%
% *var_struct.deg_in, the degree indices from the standard degrees in each 
% dimension for non square domains : 'degreeIndex'.
% 
% *var_struct.cdeg_in, the coarse degree to be used (if applicable) 
% 
% *var_struct.cheb_deg_in ,the degree indices from the standard degrees in 
% each dimension for square domains.
%
% Here the degrees can be chosen from the set [3 5 9 17 33 65 129].  
% So if 'degreeIndex', [5 5 5], the max degree of any approximate will be 
% 33 in each direction. 
    properties
        ChebRoot
        TreeGrid
        leafArray
        Errs
        Nums
        domain_in
        grid_opt = false;
        max_lengths
        GlobalMax
        LocalMax
        cheb_degs
    end
    
    methods (Abstract)
        Max = sample(obj,f,grid_opt)
        Child = splitleaf(obj,Max,set_vals)
        %[zone,domain] = splitleafGeom(zone,outerbox)
        Child = split(obj,split_dim,set_vals)
    end
    
    methods
        function obj = LSPatch(var_struct)
            
            %Call superclass constructor
            obj = obj@LeafPatch(var_struct);
            
            obj.is_geometric_refined = true;
            
            obj.max_lengths = zeros(obj.dim,1);
            obj.max_lengths(:) = inf;
            
            obj.cheb_degs = zeros(1,obj.dim);
            obj.cheb_degs(:) = 64;
            
            if isfield(var_struct,'chebDegree')
                obj.cheb_degs = var_struct.cheb_degs;
            end
            
            if isfield(var_struct,'max_lengths')
                obj.max_lengths = var_struct.max_lengths;
            end
            
            if isfield(var_struct,'domain_in')
                obj.domain_in = var_struct.domain_in;
            else
                error('An inner domain needs to be specified');
            end
            
        end
        
       %Returns structure of parameters
       function p_struct = params(obj) 
           
           p_struct.outerbox = obj.outerbox;
           p_struct.zone = obj.zone;
           p_struct.domain = obj.domain;
           p_struct.domain_in = obj.domain_in;
           p_struct.split_flag = obj.split_flag;
           p_struct.max_lengths = obj.max_lengths;
           p_struct.tol = obj.tol;
           p_struct.degs = obj.degs;
       end
       
    end
        
        
        
end