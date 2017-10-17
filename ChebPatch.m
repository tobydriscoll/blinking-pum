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
        bump
        ind_start
    end
    
    properties (Access = protected)
        deg_in %index for the standard degrees
        cdeg_in %index for the course degrees
        swap_deg_in
        split_flag %array indicating if we will split along a dimension
    end
    
    properties (Constant)
        standard_variables = load('cheb_points_matrices.mat');
        standard_degs = [3 5 9 17 33 65 129];
        invf = @(x,dom) 2/diff(dom)*x-sum(dom)/diff(dom); %takes points from a domain to [-1 1]
        forf = @(x,dom) 0.5*diff(dom)*x+0.5*sum(dom); %takes points from [-1 1] to a domain
        cheb_bump = chebfun( {0,@(x) exp(1-1./(1-x.^2)),0},[-20 -1 1 20]);
    end
    
    methods
        % Construct for the ChebPatch
        %
        %  Input:
        % domain: (dim x 2) array indiciating array for the hypercube.
        %   degs: (dim x 1) array indicating the degree of the polynomial
        %         along the dimensions.
        function obj = ChebPatch(region,zone,outerbox,deg_in,split_flag,tol,cdeg_in)
            obj.outerbox = outerbox;
            obj.region = region;
            obj.zone = zone;
            obj.domain = region;
            obj.is_geometric_refined = true; %square is always geometrically refined
            [obj.dim,~] = size(obj.domain);
            
           if nargin < 4
                obj.deg_in = zeros(1,obj.dim);
                obj.cdeg_in = zeros(1,obj.dim);
                obj.deg_in(:) = 7;
                obj.cdeg_in(:) = 3;
                obj.split_flag = ones(obj.dim,1);
                obj.tol = 1e-12;
            elseif nargin < 5
                obj.deg_in = deg_in;
                obj.cdeg_in = zeros(1,obj.dim);
                 obj.cdeg_in(:) = 3;
                obj.split_flag = ones(obj.dim,1);
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
            obj.cheb_length = prod(obj.degs);
            obj.is_leaf = true;
            obj.is_refined = false;
            obj.is_geometric_refined = true;
            
            obj.bump = cell(3,1);
            
            for k=1:obj.dim
                if obj.domain(k,1) == obj.outerbox(k,1)
                    w = @(x) obj.cheb_bump((obj.invf(x,obj.domain(k,:))+1)/2);
                elseif obj.domain(k,2) == obj.outerbox(k,2)
                    w = @(x) obj.cheb_bump((obj.invf(x,obj.domain(k,:))-1)/2);
                else
                    w = @(x) obj.cheb_bump(obj.invf(x,obj.domain(k,:)));
                end
                obj.bump{k} = w;
            end
        end
        
        
        % Returns the length of the object
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
                
                obj.values = obj.evalfGrid(grid,1,0);
                
                obj.degs = obj.cdegs;
                obj.deg_in = obj.cdeg_in;
                
                obj.cheb_length = prod(obj.cdegs);
            end
        end
        
                % Returns the length of the object
        function Refine(obj)
            if obj.iscoarse
                
                obj.iscoarse = false;
                
                obj.degs = obj.swap_degs;
                obj.deg_in = obj.swap_deg_in;
                
                grid = obj.leafGrids();
                
                obj.degs = obj.cdegs;
                obj.deg_in = obj.cdeg_in;
                
                obj.values = obj.evalfGrid(grid,1,0);
                
                obj.degs = obj.swap_degs;
                obj.deg_in = obj.swap_deg_in;
                
                obj.cheb_length = prod(obj.degs);
            end
        end
        
        % Returns the length of the object
        function ln=length(obj)
            ln = obj.cheb_length;
        end
        
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
        %    dim: dimension we are evaluating along
        %  order: order of the derivative we are calculating
        %
        % Output:
        %     ef: (length(X) x order+1) array containing the derivatves
        %         along dimension dim from 0 to order.
        function ef = evalf(obj,X,diff_dim,order)
            
            [num_pts,~] = size(X);
            
            ef = zeros(num_pts,order+1);
            
            input{1} = obj.values;
            
             for j=2:order+1
                switch obj.dim
                    case 1
                        diff_values = 2/diff(obj.domain)*obj.standard_variables.chebmatrices{obj.deg_in(1)}*obj.values;
                    otherwise
                        diff_values = (shiftdim(obj.values,diff_dim-1));
                        perm_degs = size(diff_values);
                        diff_values = reshape(diff_values,[perm_degs(1) prod(perm_degs(2:end))]);
                        diff_values = (2/diff(obj.domain(diff_dim,:)))^(j-1)*obj.standard_variables.chebmatrices{obj.deg_in(diff_dim),j-1}*diff_values;
                        diff_values = reshape(diff_values,perm_degs);
                        diff_values = shiftdim(diff_values,obj.dim-(diff_dim-1));
                end
                input{j} = diff_values;
            end
            
            for j=1:order+1
                for i=1:size(X,1)
                    G = input{j};
                    for k=1:obj.dim
                        point = obj.invf(X(i,k),obj.domain(k,:));
                        G = bary(point,G,obj.standard_variables.chebpoints{obj.deg_in(k)},obj.standard_variables.chebweights{obj.deg_in(k)});
                    end
                    ef(i,j) = G;
                end
            end
            
        end
        
        
        function ef = evalfBump(obj,X)
            
            ef = ones(size(X,1),1);
            for k=1:obj.dim
                ef = ef.*obj.bump{k}(X(:,k));
            end
            
        end
        
        function ef = evalfGridBump(obj,X)
            
            W = cell(3,1);
            
            for i=1:obj.dim
                W{i} = obj.bump{i}(X{i});
            end
            
            if obj.dim==2
                ef = W{1}*W{2}.';
            else
                ef = reshape(W{3},1,1,length(W{3})).*(W{2}'.*W{1});
            end
            
        end
        
        
        function M = interpMatrixPoints(obj,X)
            
            M = zeros(size(X,1),length(obj));
            
            for i=1:size(X,1)
                M(i,:) = obj.interpMatrixGrid(num2cell(X(i,:)));
            end
            
        end
        
        
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
        
        function ef = evalfGrid(obj,X,diff_dim,order)
            
            grid_lengths = cellfun(@(x)length(x),X);
            
            ef = zeros([grid_lengths order+1]);
            
            input{1} = obj.values;
            
            for j=2:order+1
                switch obj.dim
                    case 1
                        diff_values = 2/diff(obj.domain)*obj.standard_variables.chebmatrices{obj.deg_in(1),j-1}*obj.values;
                    otherwise
                       D = (2/diff(obj.domain(diff_dim,:)))^(j-1)*obj.standard_variables.chebmatrices{obj.deg_in(diff_dim),j-1};
                       diff_values = chebfun3t.txm(obj.values,D,diff_dim);     
                end
                input{j} = diff_values;
            end
            
            for j=1:order+1
                G = input{j};
                for k=1:obj.dim
                    %Shift the points to the right domain
                    points = obj.invf(X{k},obj.domain(k,:));
                    
                    f = bary(points,eye(obj.degs(k)),obj.standard_variables.chebpoints{obj.deg_in(k)},obj.standard_variables.chebweights{obj.deg_in(k)});
                    
                    if length(X{k})==1
                        f = f.';
                    end
                    G = chebfun3t.txm(G, f, k);
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
        function Child = splitleaf(obj,set_vals)
            
            if nargin==1
                set_vals = false;
            end
            
            lens = zeros(obj.dim,1);
            
            if obj.dim==2
                

                pref = chebfunpref();
                pref.chebfuneps = obj.tol;
                
                simple_2D_coeffs = chebfun2.vals2coeffs(obj.values);             
%                 if obj.split_flag(1)
%                     colChebtech = sum(abs(simple_2D_coeffs), 2);
%                     fCol = chebtech2({[], colChebtech});
%                     [isHappyX, cutoffX2] = happinessCheck(fCol, [], [], [], pref);
%                     lens(1) = cutoffX2+~isHappyX;
%                 end
%                 
%                 if obj.split_flag(2)
%                     rowChebtech = sum(abs(simple_2D_coeffs.'), 2);
%                     fRow = chebtech2({[], rowChebtech});
%                     [isHappyY, cutoffY2] = happinessCheck(fRow, [], [], [], pref);
%                     lens(2) = cutoffY2+~isHappyY;
%                 end
                

                if obj.split_flag(1)
                    fCol = chebtech2(obj.values);
                    [isHappyX, cutoffX2] = happinessCheck(fCol, [], [], [], pref);
                    lens(1) = cutoffX2+~isHappyX;
                end
                
                if obj.split_flag(2)
                    fRow = chebtech2(obj.values');
                    [isHappyY, cutoffY2] = happinessCheck(fRow, [], [], [], pref);
                    lens(2) = cutoffY2+~isHappyY;
                end
                
                
            else
                pref = chebfunpref();
                pref.chebfuneps = obj.tol;
                simple_3D_coeffs = chebfun3t.vals2coeffs(obj.values);
                
                colChebtech = sum(chebfun3t.unfold(abs(simple_3D_coeffs), 1), 2);
                fCol = chebtech2({[], colChebtech});
                [isHappyX, cutoffX2] = happinessCheck(fCol, [], [], [], pref);
                lens(1) = cutoffX2+~isHappyX;
                
                rowChebtech = sum(chebfun3t.unfold(abs(simple_3D_coeffs), 2), 2);
                fRow = chebtech2({[], rowChebtech});
                [isHappyY, cutoffY2] = happinessCheck(fRow, [], [], [], pref);
                lens(2) = cutoffY2+~isHappyY;
                
                
                tubeChebtech = sum(chebfun3t.unfold(abs(simple_3D_coeffs), 3), 2);
                fTube = chebtech2({[], tubeChebtech});
                [isHappyZ, cutoffZ2] = happinessCheck(fTube, [], [], [], pref);
                lens(3) = cutoffZ2+~isHappyZ;
                
            end
            
            for i=1:obj.dim
                
                %We first see if the interpolant is refined along dimension
                %i
                if obj.split_flag(i) && lens(i)<obj.degs(i)
                    obj.split_flag(i) = false;
                    
                    %If it is, find the smallest degree we can use.
                    k = find(lens(i)<=obj.standard_degs,1);
                    
                    %Slice the values for the new dimension k along i
                    if k<obj.deg_in(i)
                        obj.values = obj.slice(obj.values,1:2^(obj.deg_in(i)-k):obj.standard_degs(obj.deg_in(i)),i);
                    end
                    
                    obj.deg_in(i) = k;
                    obj.degs(i) = obj.standard_degs(k);
                    
                end
            end
            
            %We find the longest length of the box, and split along that
            %dimension
            lengths = diff(obj.domain');
            lengths(~obj.split_flag) = 0;
            [~,split_dim] = max(lengths);
            
            if ~any(obj.split_flag)
                %The leaf is refined, so return it.
                obj.is_refined = true;
                obj.cheb_length = prod(obj.degs);
                Child = obj;
            else
                
                %Return the PUPatch with the new children
                Child = obj.split(split_dim,set_vals);
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
        function Child = split(obj,split_dim,set_vals)
            
                if nargin == 2
                    set_vals = false;
                end
                
                Children = cell(1,2);
                
                %The width of the overlap
                delta = 0.5*obj.overlap*diff(obj.zone(split_dim,:));
                
                zone0 = obj.zone;
                zone1 = obj.zone;
                
                region0 = obj.region;
                region1 = obj.region;
                
                m = sum(obj.zone(split_dim,:))/2;
                
                zone0(split_dim,:) = [obj.zone(split_dim,1) m];
                zone1(split_dim,:) = [m obj.zone(split_dim,2)];
                
                region0(split_dim,:) = [max(obj.outerbox(split_dim,1),obj.zone(split_dim,1)-delta) m+delta];
                region1(split_dim,:) = [m-delta,min(obj.outerbox(split_dim,2),obj.zone(split_dim,2)+delta)];
                
                overlap_in = [m-delta,m+delta];
                
                Children{1} = ChebPatch(region0,zone0,obj.outerbox,obj.deg_in,obj.split_flag,obj.tol,obj.cdeg_in);
                
                Children{2} = ChebPatch(region1,zone1,obj.outerbox,obj.deg_in,obj.split_flag,obj.tol,obj.cdeg_in);
                
                region = obj.region;
                
                region(split_dim,:) = [region0(split_dim,1) region1(split_dim,2)];
                
                %Return the PUPatch with the new children
                Child = PUPatch(region,obj.zone,overlap_in,length(Children{1})+length(Children{2}),Children,split_dim,obj.index);
                
                if set_vals
                    for k=1:2
                        Child.children{k}.sample(obj.evalfGrid(Child.children{k}.leafGrids(),1,0));
                    end
                end
        end
        
        
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
        
        
        function plotdomain(obj)
            
            if obj.dim==2
                hold on;
                lengths = [diff(obj.domain(1,:));diff(obj.domain(2,:))];
                rectangle('position',[obj.domain(:,1)' lengths'],'LineWidth',2);
                plot(mean(obj.domain(1,:)),mean(obj.domain(2,:)),'.','MarkerSize',10,'Color','black');
                hold off;
            elseif obj.dim==3
                hold on;
                lengths = [diff(obj.domain(1,:));diff(obj.domain(2,:));diff(obj.domain(3,:))];
                center = sum(obj.domain,2)/2;
                %Vertices for Line Cube. Order matters
                X = [-1 -1 1 1 -1 -1 1 1 1 1 1 1 -1 -1 -1 -1 -1]';
                Y = [-1 1 1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1]';
                Z = [-1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 1 -1 1 1 -1]';
                %Example two cube matrix. Unit cube and one scaled/translated cube
                X1 = X*lengths(1)/2+center(1);
                Y1 = Y*lengths(2)/2+center(2);
                Z1 = Z*lengths(3)/2+center(3);
                %Single plot command for all 'cube lines'
                plot3(X1,Y1,Z1,'color','black');
                hold off;
            end
        end
        
        function plotzone(obj)
            
            if obj.dim==2
                hold on;
                lengths = [diff(obj.zone(1,:));diff(obj.zone(2,:))];
                rectangle('position',[obj.zone(:,1)' lengths'],'LineWidth',2);
                hold off;
            elseif obj.dim==3
                hold on;
                lengths = [diff(obj.zone(1,:));diff(obj.zone(2,:));diff(obj.zone(3,:))];
                center = sum(obj.zone,2)/2;
                %Vertices for Line Cube. Order matters
                X = [-1 -1 1 1 -1 -1 1 1 1 1 1 1 -1 -1 -1 -1 -1]';
                Y = [-1 1 1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1]';
                Z = [-1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 1 -1 1 1 -1]';
                %Example two cube matrix. Unit cube and one scaled/translated cube
                X1 = X*lengths(1)/2+center(1);
                Y1 = Y*lengths(2)/2+center(2);
                Z1 = Z*lengths(3)/2+center(3);
                %Single plot command for all 'cube lines'
                plot3(X1,Y1,Z1,'color','black');
                hold off;
            end
        end
        
        function IsGeometricallyRefined = IsGeometricallyRefined(obj)
            IsGeometricallyRefined = true;
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