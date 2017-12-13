classdef Chebfun2Patch<LeafPatch
    % This class handles the operations for a Chebyshev grid interpolation on a hypercube.
    
    properties
        degs %array of degrees along the dimensions
        cheb2 %grid of values to be used for interpolation
        cdegs %array of coarse degrees along the dimensions
        swap_degs %temp holder for degs
        iscoarse = false;
        Crank
        linOp
        ClinOp
        Binterp
        ind_start
    end
    
    properties (Access = protected)
        swap_deg_in
    end
    
    properties (Constant)
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
        function obj = Chebfun2Patch(domain,zone,outerbox,rank,degs,tol)

            obj = obj@LeafPatch(domain,zone,outerbox);
            
            
            if nargin < 4
                obj.Crank = 10;
                obj.degs = [64 64];
                obj.tol = 1e-12;
            elseif nargin < 5
                obj.Crank = rank;
                obj.degs = [64 64];
                obj.tol = 1e-12;
            elseif nargin < 6
                obj.Crank = rank;
                obj.degs = degs;
                obj.tol = 1e-12;
            elseif nargin < 7
                obj.Crank = rank;
                obj.degs = [64 64];
                obj.tol = tol;
            end
            
            obj.cheb_length = obj.Crank*sum(obj.degs);
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
        
        % Evaluates the approximant and its derivatives.
        %
        %  Input:
        %      X: set of points to evaluate at
        %
        % Output:
        %     ef: length(X) array containing the interpolated
        function ef = evalf(obj,X)
            
            [num_pts,~] = size(X);
            
            ef = zeros(num_pts,1);
            
            for i=1:num_pts
                ef(i) = evalfGrid(obj,num2cell(X(i,:)));
            end
            
        end
        
        % Evaluates the approximant on a grid.
        %
        %  Input:
        %      X: cellarray of grids to be evaluated on.
        %
        % Output:
        %     ef: matrix of dim(X) containing the interpolated values
        function ef = evalfGrid(obj,X)
            [X,Y] = ndgrid(X{:});
            ef = obj.cheb2(X,Y);
            
            if size(ef,1)~=size(X,1) || size(ef,2)~=size(X,2)
                ef = ef.';
            end
        end
        
        
        % Sets the values to be used for interpolation
        %
        % Input:
        %     f: values sampled at obj.points.
        function [Max] = sample(obj,f)
           % chebfunpref.setDefaults({'cheb2Prefs','maxRank'},obj.Crank);
            %Just assume we sample f for right now.
            warning('off','all');
            obj.cheb2 = chebfun2(f,[obj.domain(1,:) obj.domain(2,:)],obj.degs);
            warning('on','all');
%            Max = max(abs(obj.values(:)));
            Max = inf;
           % chebfunpref.setDefaults({'cheb2Prefs','maxRank'},512);
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
            
            r = rank(obj.cheb2);
            
            if r<obj.Crank
                
                SF = simplify(obj.cheb2);
                
                [m,n] = length(SF);
                
                if m<obj.degs(1) && n<obj.degs(2)
                
                obj.cheb2 = SF;
                
                obj.degs(1) = m;
                obj.degs(2) = n;
                obj.is_refined = true;
                obj.Crank = r;
                
                end
            end
            
            Child = obj;
            
            if ~obj.is_refined
                %Go through and split in each dimension
                for k=1:obj.dim
                    if Child.is_leaf
                        Child = split(obj,k,false);
                    else
                        Child.split(k,false);
                    end
                end
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
    
    overlap_in = [m-delta,m+delta];
    
    Children{1} = Chebfun2Patch(region0,zone0,obj.outerbox,obj.Crank,obj.degs,obj.tol);
    
    Children{2} = Chebfun2Patch(region1,zone1,obj.outerbox,obj.Crank,obj.degs,obj.tol);
    
    domain = obj.domain;
    
    domain(split_dim,:) = [region0(split_dim,1) region1(split_dim,2)];
    
    %Return the PUPatch with the new children
    Child = PUPatch(domain,obj.zone,overlap_in,length(Children{1})+length(Children{2}),Children,split_dim,obj.index);
    
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
        str = strcat(str,sprintf('rank %d ',obj.Crank));
        
        str = strcat(str,sprintf(' length %d', length(obj)));
    end
    
    
    function IsGeometricallyRefined = IsGeometricallyRefined(obj)
        IsGeometricallyRefined = true;
    end
    
    function pts = points(obj)
        pts = [];
    end
    
    
    
    function ln = length(obj)
        ln = obj.Crank*(sum(obj.degs));
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