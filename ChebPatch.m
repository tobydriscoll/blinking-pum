classdef ChebPatch<LeafPatch
    % This class handles the operations for a Chebyshev grid interpolation on a hypercube.
    
    properties
        degs %array of degrees along the dimensions
        values %grid of values to be used for interpolation
        cdegs %array of coarse degrees along the dimensions
        swap_degs %temp holder for degs
        iscoarse = false;
        linOp
        ClinOp
        Binterp
        ind_start
        deg_in %index for the standard degrees
        orig_degs
        orig_deg_in
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
        % Construct for the ChebPatch
        %
        %      Input:
        %       zone: (dim x 2) array indiciating array for the zone.
        %     domain: (dim x 2) array indiciating array for the domain.
        %   outerbox: (dim x 2) array indiciating array for the outerbox.
        %    degs_in: (dim x 1) integer array indicating the degree in each
        %                       dimension. Here the standard degrees are
        %                       [3 5 9 17 33 65 129].
        % split_flag: (dim x 1) boolean array indicating if the patch can
        %                       be split in any given dimension.
        function obj = ChebPatch(domain,zone,outerbox,deg_in,split_flag,tol,cdeg_in)
            
            %Call superclass constructor
            obj = obj@LeafPatch(domain,zone,outerbox);
            
            if nargin < 4
                obj.deg_in = zeros(1,obj.dim);
                obj.cdeg_in = zeros(1,obj.dim);
                obj.deg_in(:) = 7;
                obj.cdeg_in(:) = 3;
                obj.split_flag = true(obj.dim,1);
                obj.tol = 1e-12;
            elseif nargin < 5
                obj.deg_in = deg_in;
                obj.cdeg_in = zeros(1,obj.dim);
                obj.cdeg_in(:) = 3;
                obj.split_flag = true(obj.dim,1);
                obj.tol = 1e-12;
            elseif nargin < 6
                obj.deg_in = deg_in;
                obj.cdeg_in = zeros(1,obj.dim);
                obj.cdeg_in(:) = 3;
                obj.split_flag = split_flag;
                obj.tol = 1e-12;
            elseif nargin < 7
                obj.deg_in = deg_in;
                obj.cdeg_in = zeros(1,obj.dim);
                obj.cdeg_in(:) = 3;
                obj.split_flag = split_flag;
                obj.tol = tol;
            elseif nargin < 8
                obj.deg_in = deg_in;
                obj.cdeg_in = cdeg_in;
                obj.split_flag = split_flag;
                obj.tol = tol;
            end
            
            obj.degs = obj.standard_degs(obj.deg_in);
            obj.cdegs = obj.standard_degs(obj.cdeg_in);
            
            obj.orig_degs = obj.degs;
            obj.orig_deg_in = obj.deg_in;
            
            obj.cheb_length = prod(obj.degs);
            obj.is_refined = false;
            obj.is_geometric_refined = true;
            
        end
        
        
        % Construct for the ChebPatch
        %
        % This method coarsens the patch, resampling the patches values
        % on a coarse grid.
        function Coarsen(obj)
            if ~obj.iscoarse
                
                obj.iscoarse = true;
                
                obj.swap_degs = obj.degs;
                obj.swap_deg_in = obj.deg_in;
                
                obj.degs = obj.cdegs;
                obj.deg_in = obj.cdeg_in;
                
                grid = obj.leafGrids();
                
                obj.degs = obj.swap_degs;
                obj.deg_in = obj.swap_deg_in;
                
                obj.values = obj.evalfGrid(grid);
                
                obj.degs = obj.cdegs;
                obj.deg_in = obj.cdeg_in;
                
                obj.cheb_length = prod(obj.cdegs);
            end
        end
        
        
        % Construct for the ChebPatch
        %
        % This method refines the patch, resampling the patches values
        % on a fine grid.
        function Refine(obj)
            if obj.iscoarse
                
                obj.iscoarse = false;
                
                obj.degs = obj.swap_degs;
                obj.deg_in = obj.swap_deg_in;
                
                grid = obj.leafGrids();
                
                obj.degs = obj.cdegs;
                obj.deg_in = obj.cdeg_in;
                
                obj.values = obj.evalfGrid(grid);
                
                obj.degs = obj.swap_degs;
                obj.deg_in = obj.swap_deg_in;
                
                obj.cheb_length = prod(obj.degs);
            end
        end
        
        % Returns the length of the object
        function ln=length(obj)
            ln = obj.cheb_length;
        end
        
        % Returns the values of the object
        function vals = Getvalues(obj)
            vals = obj.values(:);
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
        
        % Returns a cell array of the grids of the domain
        function grid = leafGrids(obj)
            grid = cell(1,obj.dim);
            
            for i=1:obj.dim
                grid{i} = obj.standard_variables.chebpoints{obj.deg_in(i)};
                grid{i} = obj.forf(grid{i},obj.domain(i,:));
            end
        end
        
        % Evaluates the approximant and its derivatives.
        %
        %  Input:
        %      X: set of points to evaluate at
        %
        % Output:
        %     ef: length(X) array containing the interpolated
        function ef = evalf(obj,X,G)
            
            if nargin<3
                G = obj.values;
            end
            
            [num_pts,~] = size(X);
            
            ef = zeros(num_pts,1);
            
            for i=1:num_pts
                ef(i) = evalfGrid(obj,num2cell(X(i,:)),G);
            end
            
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
            
            G = chebfun3t.txm(obj.values, obj.standard_variables.matrices{obj.deg_in(diff_dim)}, order)/diff(obj.domain(diff_dim,:))^order;
            
            if nargin<4
                ef = G(:);
            else
                ef = evalfDiffGrid(obj,diff_dim,order,X,G);
            end
        end
        
        % Evaluates the approximant on a grid.
        %
        %  Input:
        %      X: cellarray of grids to be evaluated on.
        %
        % Output:
        %     ef: matrix of dim(X) containing the interpolated values
        function ef = evalfGrid(obj,X,G)
            
            if nargin<3
                G = obj.values;
            end
            
            for k=1:obj.dim
                %Shift the points to the right domain
                points = obj.invf(X{k},obj.domain(k,:));
                
                f = bary(points,eye(obj.degs(k)),obj.standard_variables.chebpoints{obj.deg_in(k)},obj.standard_variables.chebweights{obj.deg_in(k)});
                
                if length(X{k})==1
                    f = f.';
                end
                
                G = chebfun3t.txm(G, f, k);
            end
            ef = G;
        end
        
        % Evaluates the derivative on a grid.
        %
        %  Input:
        %      diff_dim: dimension derivative is taken in
        %         order: order of derivative (up to 2 right now)
        %             X: cellarray of grids to be evaluated on.
        %
        % Output:
        %     ef: matrix of dim(X) containing the interpolated derivative values
        function ef = evalfDiffGrid(obj,diff_dim,order,X)
            if nargin<3
                order = 1;
            end
            
            G = chebfun3t.txm(obj.values, obj.standard_variables.chebmatrices{obj.deg_in(diff_dim),order}, diff_dim)/diff(obj.domain(diff_dim,:))^order;
            
            if nargin<4
                ef = G;
            else
                ef = evalfGrid(obj,X,G);
            end
        end
        
        %  interpMatrixPoints(obj,X)
        %  This method creates a interpolating matrix given a list of
        %  points.
        %
        %  Input:
        %      X: list of points.
        %
        % Output:
        %      M: interpolating matrix.
        function M = interpMatrixPoints(obj,X)
            
            M = zeros(size(X,1),length(obj));
            
            for i=1:size(X,1)
                M(i,:) = obj.interpMatrixGrid(num2cell(X(i,:)));
            end
            
        end
        
        %  interpMatrixGrid(obj,grid)
        %  This method creates a interpolating matrix given a grid.
        %
        %  Input:
        %      X: cell array of grid values.
        %
        % Output:
        %      M: interpolating matrix.
        function M = interpMatrixGrid(obj,grid)
            G = obj.leafGrids();
            if obj.dim==1
                if iscell(grid)
                    M = barymat(grid{1},G{1});
                else
                    M = barymat(grid,G{1});
                end
            elseif obj.dim==2
                M = kron(barymat(grid{2},G{2}),barymat(grid{1},G{1}));
            elseif obj.dim==3
                M = kron(barymat(grid{3},G{3}),kron(barymat(grid{2},G{2}),barymat(grid{1},G{1})));
            end
        end
        
        % Sets the values to be used for interpolation
        %
        % Input:
        %     f: values sampled at obj.points.
        function [Max] = sample(obj,f,grid_opt)
            
            if(nargin==2)
                grid_opt = false;
            end
            
            %Just assume we sample f for right now.
            if ~isnumeric(f)
                
                if ~grid_opt
                    if obj.dim==1
                        obj.values = f(obj.points());
                    else
                        obj.values = reshape(f(obj.points()),obj.degs);
                    end
                else
                    %If a function is more efficient on a grid, evaluate it
                    %as so.
                    obj.values = f(obj.leafGrids());
                end
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
            Max = max(abs(obj.values(:)));
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
        function Child = splitleaf(obj,Max,set_vals)
            
            if nargin==2
                set_vals = false;
            end
            
            vscale = Max;
            
            loc_tol = obj.tol^(7/8);
            
            cutoff = zeros(obj.dim,1);
            
            local_max = max(abs(obj.values(:)));
            
            if obj.dim==1
                fCol = chebtech2(obj.values);
                hscale = diff(obj.domain);
                tol = loc_tol*max(vscale./local_max,hscale);
                cutoff(1) = length(simplify(fCol, tol))+1;
            else
                
                if obj.dim==2
                    coeffs = chebfun2.vals2coeffs(obj.values);
                else
                    coeffs = chebfun3.vals2coeffs(obj.values);
                end
                
                for k=1:obj.dim
                    
                    if obj.split_flag(k)
                        
                        colChebtech = chebfun3t.unfold(coeffs, k);
                        colChebtech = sum(abs(colChebtech),2);
                        fCol = chebtech2({[],colChebtech});
                        hscale = diff(obj.domain(k,:));
                        
                        tol = loc_tol*max(vscale./local_max,hscale);
                        cutoff(k) = length(simplify(fCol, tol))+1;
                        
                    end
                    
                end
            end
            
            for k=1:obj.dim
                sliceSample(obj,k,cutoff(k));
            end
            
            
            
            
            if ~any(obj.split_flag)
                %The leaf is refined, so return it.
                obj.is_refined = true;
                obj.cheb_length = prod(obj.degs);
                Child = obj;
            else
                
%                 ind = 1:obj.dim;
%                 ind = ind(obj.split_flag);
%                 
%                 [~,split_dim] = max(diff(obj.domain(obj.split_flag,:).',1));
%                 split_dim = ind(split_dim);
%                 
%                 Child = split(obj,split_dim,set_vals);
                
                
                Child = obj;
                %Go through and split in each unresolved direction
                for k=1:obj.dim
                    if obj.split_flag(k)
                        if Child.is_leaf
                            Child = split(obj,k,set_vals);
                        else
                            Child.split(k,set_vals);
                        end
                    end
                end
                
                
                
            end
        end
        
        % The method will slice a sample in a given dimension.
        %
        %     Input:
        %      d_in: splitting dimension
        %       len: length to be sliced too
        function [] = sliceSample(obj,d_in,len)
            if obj.split_flag(d_in) && len<obj.degs(d_in)
                obj.split_flag(d_in) = false;
                
                %If it is, find the smallest degree we can use.
                k = find(len<=obj.standard_degs,1);
                
                if obj.dim==1
                    if k<obj.deg_in(d_in)
                        obj.values = obj.values(1:2^(obj.deg_in(d_in)-k):obj.standard_degs(obj.deg_in(d_in)));
                    end
                else
                    
                    %Slice the values for the new dimension k along i
                    if k<obj.deg_in(d_in)
                        obj.values = obj.slice(obj.values,1:2^(obj.deg_in(d_in)-k):obj.standard_degs(obj.deg_in(d_in)),d_in);
                    end
                end
                obj.deg_in(d_in) = k;
                obj.degs(d_in) = obj.standard_degs(k);
                
                obj.cheb_length = prod(obj.degs);
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
            
            Children = cell(1,2);
            
            %The width of the overlap
            delta = 0.5*obj.overlap*diff(obj.zone(split_dim,:));
            
            zone0 = obj.zone;
            zone1 = obj.zone;
            
            region0 = obj.domain;
            region1 = obj.domain;
            
            m = sum(obj.zone(split_dim,:))/2;
            
            zone0(split_dim,:) = [obj.zone(split_dim,1) m];
            zone1(split_dim,:) = [m obj.zone(split_dim,2)];
            
            region0(split_dim,:) = [max(obj.outerbox(split_dim,1),obj.zone(split_dim,1)-delta) m+delta];
            region1(split_dim,:) = [m-delta,min(obj.outerbox(split_dim,2),obj.zone(split_dim,2)+delta)];
            
            
            Children{1} = ChebPatch(region0,zone0,obj.outerbox,obj.deg_in,obj.split_flag,obj.tol,obj.cdeg_in);
            Children{2} = ChebPatch(region1,zone1,obj.outerbox,obj.deg_in,obj.split_flag,obj.tol,obj.cdeg_in);
            
            Children{1}.orig_degs = obj.orig_degs;
            Children{1}.orig_deg_in = obj.orig_deg_in;
            
            Children{2}.orig_degs = obj.orig_degs;
            Children{2}.orig_deg_in = obj.orig_deg_in;
            
            domain = obj.domain;
            
            domain(split_dim,:) = [region0(split_dim,1) region1(split_dim,2)];
            
            %Return the PUPatch with the new children
            Child = PUPatch(domain,obj.zone,length(Children{1})+length(Children{2}),Children,split_dim,obj.index);
            
            if set_vals
                for k=1:2
                    Child.children{k}.sample(obj.evalfGrid(Child.children{k}.leafGrids()));
                end
            end
        end
        
        % String function for object.
        %
        %     Input:
        %
        %     Child: string of object.
        function str = toString(obj)
            str = '';
            for i=1:obj.dim
                str = strcat(str,sprintf(' [%0.3f %0.3f]',obj.domain(i,1),obj.domain(i,2)));
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
        
        function IsGeometricallyRefined = IsGeometricallyRefined(obj)
            IsGeometricallyRefined = true;
        end
        
        
        
        function reset(obj)
            obj.is_refined = false;
            obj.is_geometric_refined = false;
            obj.degs = obj.orig_degs;
            obj.deg_in = obj.orig_deg_in;
            obj.split_flag = ones(obj.dim,1);
        end
        
    end
    
    methods (Static)
        %This method slices an array along a certain dimension
        function out = slice(A, ix, dim)
            subses = repmat({':'}, [1 ndims(A)]);
            subses{dim} = ix;
            out = A(subses{:});
        end
    end
    
end