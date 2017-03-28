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
        standard_variables = load('cheb_points_matrices.mat');
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
        function ef = evalf(obj,X,diff_dim,order)
            
            ef = zeros(length(X),order+1);
            
            input{1} = obj.values;
            
            if order>0
                switch obj.dim
                    case 1
                        diff_values = 2*(obj.domain(2)-obj.domain(1))^-1*obj.standard_variables.chebmatrices{obj.deg_in(1)}*obj.values;
                    otherwise
                        diff_values = shiftdim(obj.values,diff_dim-1);
                        perm_degs = size(diff_values);
                        
                        if obj.dim>2
                            diff_values = reshape(diff_values,[perm_degs(1) prod(perm_degs(2:end))]);
                        end
                        
                        diff_values = 2*(obj.domain(diff_dim,2)-obj.domain(diff_dim,1))^-1*...
                            obj.standard_variables.chebmatrices{obj.deg_in(diff_dim)}*diff_values;
                        
                        if obj.dim>2
                            diff_values = reshape(diff_values,perm_degs);
                        end
                        
                        diff_values = shiftdim(diff_values,obj.dim-(diff_dim-1));
                end
               
               input{2} = diff_values;
           end
            
            for j=1:order+1
                for i=1:size(X,1)
                    G = input{j};
                    for k=1:obj.dim
                        
                        point = 2/(obj.domain(k,2)-obj.domain(k,1))*X(i,k)-(obj.domain(k,2)+obj.domain(k,1))/(obj.domain(k,2)-obj.domain(k,1));
                        G = bary(point,G,obj.standard_variables.chebpoints{obj.deg_in(k)},obj.standard_variables.chebweights{obj.deg_in(k)});
                        
                    end
                    ef(i,j) = G;
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
        
        function fx = bary(x, fvals, deg_ind)
            %BARY   Barycentric interpolation formula.
            %   BARY(X, FVALS, XK, VK) uses the 2nd form barycentric formula with weights VK
            %   to evaluate an interpolant of the data {XK, FVALS(:,k)} at the points X.
            %   Note that XK and VK should be column vectors, and FVALS, XK, and VK should
            %   have the same length.
            %
            %   BARY(X, FVALS) assumes XK are the 2nd-kind Chebyshev points and VK are the
            %   corresponding barycentric weights.
            %
            %   If size(FVALS, 2) > 1 then BARY(X, FVALS) returns values in the form
            %   [F_1(X), F_2(X), ...], where size(F_k(X)) = size(X) for each k.
            %
            %
            % See also CHEBTECH.CLENSHAW.
            
            % Copyright 2016 by The University of Oxford and The Chebfun Developers.
            % See http://www.chebfun.org/ for Chebfun information.
            
            % Parse inputs:
            
            lengths = size(fvals);
            
            n = lengths(1);
            m = prod(lengths(2:end));
            
            numdims = size(lengths,2);
            
            fvals = reshape(fvals,[n m]);
            
            sizex = size(x);
            
            xk = ChebPatch.standard_variables.chebpoints{deg_ind};
            vk = ChebPatch.standard_variables.chebweights{deg_ind};
            
            if ( ~all(sizex) )
                fx = x;
                return
            end
            
            % Check that input is a column vector:
            if ( (sizex(2) > 1) )
                x = x(:);
            end
            
            % % The function is a constant.
            % if ( n == 1 )
            %     fx = repmat(fvals, length(x), 1);
            %     return
            % end
            
            % The function is NaN.
            if ( any(isnan(fvals)) )
                fx = NaN(length(x), m);
                return
            end
            
            % The main loop:
            if ( numel(x) < 4*length(xk) )  % Loop over evaluation points
                % Note: The value "4" here was detemined experimentally.
                
                % Initialise return value:
                fx = zeros(size(x, 1), m);
                
                % Loop:
                for j = 1:numel(x),
                    xx = vk ./ (x(j) - xk);
                    fx(j,:) = (xx.'*fvals) / sum(xx);
                end
            else                            % Loop over barycentric nodes
                % Initialise:
                num = zeros(size(x, 1), m);
                denom = num;
                
                % Loop:
                for j = 1:length(xk),
                    tmp = (vk(j) ./ (x - xk(j)));
                    num = num + bsxfun(@times, tmp, fvals(j,:));
                    denom = bsxfun(@plus, denom, tmp);
                end
                fx = num ./ denom;
            end
            
            % Try to clean up NaNs:
            for k = find(isnan(fx(:,1)))'       % (Transpose as Matlab loops over columns)
                index = find(x(k) == xk, 1);    % Find the corresponding node
                if ( ~isempty(index) )
                    fx(k,:) = fvals(index,:);   % Set entry/entries to the correct value
                end
            end
            
            % Reshape the output:
            if (numel(x)==1) && numdims>2
                fx = reshape(fx, lengths(2:end));
            elseif numdims>2
                fx = reshape(fx, [numel(x) lengths(2:end)]);
            elseif numdims==2
                fx = transpose(fx);
            end
            
        end
    end
    
end