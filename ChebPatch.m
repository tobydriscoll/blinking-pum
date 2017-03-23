classdef ChebPatch<LeafPatch
    % This class handles the operations for a Chebyshev grid interpolation on a hypercube.
    
    properties
        degs %array of degrees along the dimensions
        values %grid of values to be used for interpolation
    end
    
    properties (Access = protected)
        deg_in %index for the standard degrees
        split_flag %array indicating if we will split along a dimension
    end
    
    properties (Constant)
        standard_points = load('cheb_points_matrices.mat','chebpoints');
        standard_matrices = load('cheb_points_matrices.mat','chebmatrices');
        standard_degs = [3 5 9 17 33 65];
    end
    
    methods
        % Construct for the ChebPatch
        %
        %  Input:
        % domain: (dim x 2) array indiciating array for the hypercube.
        %   degs: (dim x 1) array indicating the degree of the polynomial
        %         along the dimensions.
        function obj = ChebPatch(domain,deg_in,split_flag)
            
            obj.domain = domain;
            [obj.dim,~] = size(obj.domain);
            
            if nargin < 2
                obj.deg_in = zeros(1,obj.dim);
                obj.deg_in(:) = 6;
                obj.split_flag = ones(obj.dim,1);
            elseif nargin < 3
                obj.deg_in = deg_in;
                obj.split_flag = ones(obj.dim,1);
            else
                obj.deg_in = deg_in;
                obj.split_flag = split_flag;
            end
            
            %To do: deal with boundary
            obj.degs = obj.standard_degs(obj.deg_in);
            obj.cheb_length = prod(obj.degs);
            obj.is_leaf = true;
            obj.is_refined = false;
        end
        
        
        % Returns the length of the object
        function ln=length(obj)
            ln = obj.cheb_length;
        end
        
        % Returns the points of the function
        function pts = points(obj)
            
            C = cell(obj.dim,1);
            
            for i=1:obj.dim
                C{i} = chebpts(obj.degs(i),obj.domain(i,:));
            end
            [out{1:obj.dim}] = ndgrid(C{:});
            
            pts = zeros(numel(out{1}),obj.dim);
            
            for i=1:obj.dim
                pts(:,i) = out{i}(:);
            end
            
        end
        
        % Evaluates the approximant and its derivatives.
        %
        %  Input:
        %      X: set of points to evaluate at
        %    dim: dimension we are evaluating along
        %  order: order of the derivative we are calculating
        %
        % Output:
        %     ef: (length(X) x order+1) array containing the derivatves
        %         along dimension dim from 0 to order.
        function ef = evalf(obj,X,dim,order)
            
            
            
            ef = zeros(length(X),order+1);
            
            %Assume we are dim 2 and up
            StartG = chebfun(obj.values,obj.domain(1,:));
            
            % This is for a list of points. If we have a grid
            % we can do this more efficiently
            for i=1:length(X)
                for ord=0:order
                    
                    for k=1:obj.dim
                        
                        switch k
                            case 1
                                ChebG = StartG;
                            case obj.dim
                                ChebG = chebfun(G',obj.domain(k,:));
                            otherwise
                                ChebG = chebfun(G,obj.domain(k,:));
                        end
                        
                        
                        switch k
                            case obj.dim
                                if k==dim
                                    ef(i,ord+1) = feval(diff(ChebG,ord),X(i,k));
                                else
                                    ef(i,ord+1) = feval(ChebG,X(i,k));
                                end
                            case obj.dim-1
                                if k==dim
                                    G = feval(diff(ChebG,ord),X(i,k));
                                else
                                    G = feval(ChebG,X(i,k));
                                end
                            otherwise
                                if k==dim
                                    G = reshape(feval(diff(ChebG,ord),X(i,k)),obj.degs(k+1:end));
                                else
                                    G = reshape(feval(ChebG,X(i,k)),obj.degs(k+1:end));
                                end
                        end
                    end
                end
            end
            
        end
        
        % Sets the values to be used for interpolation
        %
        % Input:
        %     f: values sampled at obj.points.
        function sample(obj,f)
            %Just assume we sample f for right now.
            switch obj.dim
                case 1
                    [~,n2] = size(f);
                    if n2 == 1
                        obj.values = f';
                    else
                        obj.values = f;
                    end
                    
                otherwise
                    obj.values = reshape(f,obj.degs);
            end
        end
        
        % The method determines if a splitting is needed, and creates
        % the new children if splitting is required.
        %
        %     Input:
        %   overlap: overlap intended to be used for the splitting
        %
        %    Output:
        %     Child: obj if no splitting is needed. If a splitting
        %            is needed, Child is the PUPatch object with
        %            the new children.
        function Child = splitleaf(obj)
            
            G = zeros(obj.dim,1);
            
            for i=1:obj.dim
                if obj.split_flag(i)
                    %Shift F so i is the diminsion we are going along.
                    perm = circshift((1:obj.dim)',-i+1);
                    Fi = permute(obj.values,perm);
                    %Figure out deg along dim i,
                    len = length(simplify(chebfun(Fi,obj.domain(i,:))));
                    
                    G(i) = len;
                    
                end
            end
               
            for i=1:obj.dim
                    if obj.split_flag(i) && G(i)<obj.degs(i)-1
                        obj.split_flag(i) = false;
                        k = find(G(i)<=obj.standard_degs,1);
                        
                        if k<obj.deg_in(i)
                            obj.values = obj.slice(obj.values,1:2*(obj.deg_in(i)-k):obj.standard_degs(obj.deg_in(i)),i);
                        end
                        
                        obj.deg_in(i) = k;
                    end
            end
            
           
            if ~any(obj.split_flag)
                %The leaf is refined, so return it.
                obj.is_refined = true;
                Child = obj;
            else
                [~,split_dim] = max(G);
                
                children = cell(1,2);
                
                delta = 0.5*(1+PUWeights.overlap)*...
                    (obj.domain(split_dim,2)-obj.domain(split_dim,1));
                
                domain0 = obj.domain;
                domain1 = obj.domain;
                
                domain0(split_dim,:) = [obj.domain(split_dim,1) obj.domain(split_dim,1)+delta];
                domain1(split_dim,:) = [obj.domain(split_dim,2)-delta obj.domain(split_dim,2)];
                
                overlap_in = [obj.domain(split_dim,2)-delta obj.domain(split_dim,1)+delta];
                children{1} = ChebPatch(domain0,obj.deg_in,obj.split_flag);
                children{2} = ChebPatch(domain1,obj.deg_in,obj.split_flag);
                
                %Return the PUPatch with the new children
                Child = PUPatch(obj.domain,overlap_in,length(children{1})+length(children{2}),children,split_dim);
                
            end
            
        end
        
    end
    
    methods (Static)
        %This method slices an array along a certain dimension 
        %(Matlab should have a function that does this)
        function out = slice(A, ix, dim)
            subses = repmat({':'}, [1 ndims(A)]);
            subses{dim} = ix;
            out = A(subses{:});
        end
        
    end
    
end