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
        deg_in
        cheb_deg_in
        cdeg_in
        cdegs
        degs
        domain_in
        grid_opt = false;
        coeffs
        max_lengths
        GlobalMax
        LocalMax
    end
    
    properties (Constant)
        standard_variables = load('cheb_points_matrices.mat');
        standard_degs = [3 5 9 17 33 65 129];
        invf = @(x,dom) 2/diff(dom)*x-sum(dom)/diff(dom); %takes points from a domain to [-1 1]
        forf = @(x,dom) 0.5*diff(dom)*x+0.5*sum(dom); %takes points from [-1 1] to a domain
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
            
            obj.tol = 1e-7;
            obj.deg_in = zeros(obj.dim,1);
            obj.deg_in(:) = 3;
            
            obj.cheb_deg_in = zeros(obj.dim,1);
            obj.cheb_deg_in(:) = 4;
            
            obj.cdeg_in = zeros(obj.dim,1);
            obj.cdeg_in(:) = 3;
            obj.split_flag = true(obj.dim,1);
            obj.max_lengths = inf(obj.dim,1);
            
            obj.max_lengths = zeros(obj.dim,1);
            obj.max_lengths(:) = inf;
            
            
            if isfield(var_struct, 'deg_in')
                obj.deg_in = var_struct.deg_in;
            end
            
            if isfield(var_struct, 'cheb_deg_in')
                obj.deg_in = var_struct.deg_in;
            end
            
            if isfield(var_struct,'split_flag')
                obj.split_flag = var_struct.split_flag;
            end
            
            if isfield(var_struct,'max_lengths')
                obj.max_lengths = var_struct.max_lengths;
            end
            
            if isfield(var_struct,'tol')
                obj.tol = var_struct.tol;
            end
            
            if isfield(var_struct,'cdeg_in')
                obj.cdeg_in = var_struct.cdeg_in;
            end
            
            if isfield(var_struct,'domain_in')
                obj.domain_in = var_struct.domain_in;
            else
                error('An inner domain needs to be specified');
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
           p_struct.domain_in = obj.domain_in;
           p_struct.split_flag = obj.split_flag;
           p_struct.max_lengths = obj.max_lengths;
           p_struct.tol = obj.tol;
           p_struct.cdeg_in = obj.cdeg_in;
       end
        
        % TODO. Figure out what to do here!
        function ln=length(obj)
            ln = prod(obj.degs);
        end
        
        % Returns the points of the function
        function pts = points(obj)
            
            C = cell(obj.dim,1);
            
            for i=1:obj.dim
                C{i} = obj.standard_variables.chebpoints{obj.deg_in(i)};
                C{i} = obj.forf(C{i},obj.domain(i,:));
            end
            [out{1:obj.dim}] = ndgrid(C{:});
            
            pts = zeros(numel(out{1}),obj.dim);
            
            for i=1:obj.dim
                pts(:,i) = out{i}(:);
            end
            
        end
        
        % TODO. Figure out what to do here!
        function ef = evalf(obj,X,G)
            if nargin<3
                G = obj.coeffs;
            end
            
            [num_pts,~] = size(X);
            
            ef = zeros(num_pts,1);
            
            for i=1:num_pts
                ef(i) = evalfGrid(obj,num2cell(X(i,:)),G);
            end
            
        end
        
        function ef = evalfGrid(obj,X,G)
            
            if nargin<3
                G = obj.coeffs;
            end
            
            for k=1:obj.dim
                %Shift the points to the domain [-1 1]x[-1 1]
                X{k} = obj.invf(X{k},obj.domain(k,:));
                
                %Evaluate the points at the Chebyshev polynomials
                F = clenshaw(X{k},eye(obj.degs(k)));
                
                %Multiply the coefficients with F
                G = chebfun3t.txm(G, F, k);
            end
            
            ef  = G;
            
        end
        
        % Evaluates the approximant and its derivatives.
        %
        %  Input:
        %      X: set of points to evaluate at
        %
        % Output:
        %     ef: length(X) array containing the interpolated
        function ef = Diff(obj,diff_dim,order,X)
            if nargin<3
                order = 1;
            end
            
            unContractedModes = [1:diff_dim-1, diff_dim+1:obj.dim];  
            G = obj.computeDiffCoeffs(chebfun3t.unfold(obj.coeffs,diff_dim))/diff(obj.domain(diff_dim,:))^order;
            G = chebfun3t.fold(G,obj.degs,diff_dim,unContractedModes);
            
            if nargin<4
                ef = G;
            else
                ef = evalfDiffGrid(obj,diff_dim,order,X,G);
            end
        end
        
        
        
        
    end
    
    
    methods (Static)
        function cout = computeDiffCoeffs(c)
            %COMPUTEDERCOEFFS   Recurrence relation for coefficients of derivative.
            %   C is the matrix of Chebyshev coefficients of a (possibly array-valued)
            %   CHEBTECH object.  COUT is the matrix of coefficients for a CHEBTECH object
            %   whose columns are the derivatives of those of the original.
            
            [n, m] = size(c);
            cout = zeros(n, m);                        % Initialize vector {c_r}
            w = repmat(2*(1:n-1)', 1, m);
            v = w.*c(2:end,:);                           % Temporal vector
            cout(n-1:-2:1,:) = cumsum(v(n-1:-2:1,:), 1); % Compute c_{n-2}, c_{n-4}, ...
            cout(n-2:-2:1,:) = cumsum(v(n-2:-2:1,:), 1); % Compute c_{n-3}, c_{n-5}, ...
            cout(1,:) = .5*cout(1,:);                    % Adjust the value for c_0
        end
    end
end