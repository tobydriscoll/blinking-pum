classdef ChebPatch<LeafPatch
    % This class handles the operations for a Chebyshev grid interpolation on a hypercube.
    
    properties
        degs %array of degrees along the dimensions
        values %grid of values to be used for interpolation
    end
    
    properties (Access = protected)
        max_in %max index for the standard degrees
        split_flag %array indicating if we will split along a dimension
        max_deg %max degree to be used
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
        function obj = ChebPatch(domain,degs,max_in)
            
            if nargin == 2
                obj.max_in = 6;
                obj.max_deg = 65;
            else
                obj.max_in = max_in;
                obj.max_deg = obj.standard_degs(max_in);
            end
            
            %To do: deal with boundary
            obj.domain = domain;
            [obj.dim,~] = size(obj.domain);
            obj.degs = degs;
            obj.cheb_length = prod(degs);
            obj.split_flag = ones(obj.dim,1);
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
            
            % This is for a list of points. If we have a grid
            % we can do this more efficiently
            for i=1:length(X)
                for ord=0:order
                    G = obj.values;
                    for k=1:obj.dim
                        switch k
                            case obj.dim
                                if k==dim
                                    ef(i,ord+1) = feval(diff(chebfun(G',obj.domain(k,:)),ord),X(i,k));
                                else
                                    ef(i,ord+1) = feval(chebfun(G',obj.domain(k,:)),X(i,k));
                                end
                            case obj.dim-1
                                if k==dim
                                    G = feval(diff(chebfun(G,obj.domain(k,:)),ord),X(i,k));
                                else
                                    G = feval(chebfun(G,obj.domain(k,:)),X(i,k));
                                end
                            otherwise
                                if k==dim
                                    G = reshape(feval(diff(chebfun(G,obj.domain(k,:)),ord),X(i,k)),obj.degs(k+1:end));
                                else
                                    G = reshape(feval(chebfun(G,obj.domain(k,:)),X(i,k)),obj.degs(k+1:end));
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
        function Child = splitleaf(obj,overlap)
            
            G = zeros(obj.dim,1);
            
            for i=1:obj.dim
                if obj.split_flag(i)
                    %Shift F so i is the diminsion we are going along.
                    perm = circshift((1:obj.dim)',-i+1);
                    Fi = permute(obj.values,perm);
                    %Figure out deg along dim i,
                    G(i) = length(simplify(chebfun(Fi,obj.domain(i,:))));
                    
                    obj.split_flag(i) = G(i)>=obj.max_deg-1;
                end
            end
            
            
            %Figure out which dimensions are resolved.
            for i=1:obj.dim
                if ~obj.split_flag(i)
                    k = find(G(i)<=obj.standard_degs,1);
                    obj.degs(i) = obj.standard_degs(k);
                    %Set the dim to the appropriate number
                    obj.values = obj.slice(obj.values,1:2*(obj.max_in-k):obj.standard_degs(obj.max_in),i);
                end
            end
            
            if ~any(obj.split_flag)
                %The leaf is refined, so return it.
                obj.is_refined = true;
                Child = obj;
            else
                [~,split_dim] = max(G);
                
                children = cell(1,2);
                
                delta = 0.5*(1+overlap)*...
                    (obj.domain(split_dim,2)-obj.domain(split_dim,1));
                
                domain0 = obj.domain;
                domain1 = obj.domain;
                
                domain0(split_dim,:) = [obj.domain(split_dim,1) obj.domain(split_dim,1)+delta];
                domain1(split_dim,:) = [obj.domain(split_dim,2)-delta obj.domain(split_dim,2)];
                
                children{1} = ChebPatch(domain0,obj.degs);
                children{2} = ChebPatch(domain1,obj.degs);
                
                %Return the PUPatch with the new children
                Child = PUPatch(obj.domain,length(children{1})+length(children{2}),children,overlap,split_dim);
                
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