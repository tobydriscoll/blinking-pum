classdef ChebPatch<LeafPatch
    % This class handles the operations for a Chebyshev grid interpolation on a hypercube.
    
    properties
        degs %array of degrees along the dimensions
        values %grid of values to be used for interpolation
        tol %tolerence used
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
            
            obj.tol = 1e-14;
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
                %C{i} = chebpts(obj.degs(i),obj.domain(i,:));
                C{i} = obj.standard_variables.chebpoints{obj.deg_in(i)};
                C{i} = 0.5*diff(obj.domain(i,:))*C{i}+0.5*sum(obj.domain(i,:));
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
                        diff_values = (shiftdim(obj.values,diff_dim-1));
                        perm_degs = size(diff_values);
                        diff_values = reshape(diff_values,[perm_degs(1) prod(perm_degs(2:end))]);
                        diff_values = 2*(obj.domain(diff_dim,2)-obj.domain(diff_dim,1))^-1*obj.standard_variables.chebmatrices{obj.deg_in(diff_dim),1}*diff_values;
                        diff_values = reshape(diff_values,perm_degs);
                        diff_values = shiftdim(diff_values,obj.dim-(diff_dim-1));
                end
                
                input{2} = diff_values;
            end
            
            for j=1:order+1
                for i=1:size(X,1)
                    G = input{j};
                    for k=1:obj.dim
                        
                        point = 2/(diff(obj.domain(k,:)))*X(i,k)-(obj.domain(k,2)+obj.domain(k,1))/(obj.domain(k,2)-obj.domain(k,1));
                        G = bary(point,G,obj.standard_variables.chebpoints{obj.deg_in(k)},obj.standard_variables.chebweights{obj.deg_in(k)});
                        
                    end
                    ef(i,j) = G;
                end
            end
            
        end
        
        function ef = evalfGrid(obj,X,diff_dim,order)
            
            grid_lengths = cellfun(@(x)length(x),X);
            
            ef = zeros([grid_lengths order]);
            
            input{1} = obj.values;
            
            if order>0
                switch obj.dim
                    case 1
                        diff_values = 2*(obj.domain(2)-obj.domain(1))^-1*obj.standard_variables.chebmatrices{obj.deg_in(1)}*obj.values;
                    otherwise
                        diff_values = (shiftdim(obj.values,diff_dim-1));
                        perm_degs = size(diff_values);
                        diff_values = reshape(diff_values,[perm_degs(1) prod(perm_degs(2:end))]);
                        diff_values = 2*(obj.domain(diff_dim,2)-obj.domain(diff_dim,1))^-1*obj.standard_variables.chebmatrices{obj.deg_in(diff_dim),1}*diff_values;
                        diff_values = reshape(diff_values,perm_degs);
                        diff_values = shiftdim(diff_values,obj.dim-(diff_dim-1));
                end
                
                input{2} = diff_values;
            end
            
            for j=1:order+1
                G = input{j};
                for k=1:obj.dim
                    %Shift the points to the right domain
                    points = 2/(diff(obj.domain(k,:)))*X{k}-sum(obj.domain(k,:))./diff(obj.domain(k,:));
                    
                    G = bary(points,G,obj.standard_variables.chebpoints{obj.deg_in(k)},obj.standard_variables.chebweights{obj.deg_in(k)});
                    G = shiftdim(G,1);
                end
                
                if order>0
                    subses = repmat({':'}, [1 ndims(ef)]);
                    subses{obj.dim+1} = j;
                    ef(subses{:}) = G;
                else
                    ef = G;
                end
            end
            
            
        end
        
        % Sets the values to be used for interpolation
        %
        % Input:
        %     f: values sampled at obj.points.
        function sample(obj,f)
            %Just assume we sample f for right now.
            if ~isnumeric(f)
                obj.values = reshape(f(obj.points()),obj.degs);
            else
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
                %if obj.split_flag(i)
                Fi = shiftdim(obj.values,i-1);
                sizes = size(Fi);
                Fi = reshape(Fi,[sizes(1) prod(sizes(2:end))]);
                %Figure out deg along dim i,
                len = simplify(Fi,obj.tol);
                G(i) = len;
                %end
            end
            
            max_in = 0;
            split_dim = 1;
            for i=1:obj.dim
                
                %                 if diff(obj.domain(i,:))>max_in
                %                     max_in = diff(obj.domain(i,:));
                %                     split_dim = i;
                %                 end
                
                if G(i)>max(max_in,obj.degs(i)-1) && obj.split_flag(i)
                    max_in = G(i);
                    split_dim = i;
                end
                
                if obj.split_flag(i) && G(i)<obj.degs(i)
                    obj.split_flag(i) = false;
                    k = find(G(i)<=obj.standard_degs,1);
                    
                    if k<obj.deg_in(i)
                        obj.values = obj.slice(obj.values,1:2*(obj.deg_in(i)-k):obj.standard_degs(obj.deg_in(i)),i);
                    end
                    
                    obj.deg_in(i) = k;
                    obj.degs(i) = obj.standard_degs(k);
                    
                end
            end
            
            
            if ~any(obj.split_flag)
                %The leaf is refined, so return it.
                obj.is_refined = true;
                obj.cheb_length = prod(obj.degs);
                Child = obj;
            else
                
                %[~,split_dim] = max(G);
                
                
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
        
        function str = toString(obj)
            str = '';
            for i=1:obj.dim
                str = strcat(str,sprintf(' [%0.3f %0.3f]',obj.domain(i,1),obj.domain(i,1)));
                if i<obj.dim
                    str = strcat(str,' x ');
                end
            end
            str = strcat(str,' degrees: ');
            for i=1:obj.dim
                str = strcat(str,sprintf(' %d ',obj.degs(i)));
            end
            
            str = strcat(str,sprintf(' length %d', obj.length));
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