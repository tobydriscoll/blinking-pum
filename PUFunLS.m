classdef PUFunLS < handle & matlab.mixin.Copyable
    % PUFun2DLS PUFun class for representing n-d functions on non-square domains.
    %
    % This class represents smooth multivariate functions on non-square domains
    % with a partition of unity approximation. This class
    % automatically finds a set of overlapping domains that are adapted to the
    % features of the function, and the blends locally defined Chebyshev
    % approximations on the domains with a partition of unity.
    %
    % PUFun2DLS(f,domain_in,domain_out) constructs a partition of unity
    % approximation representing f on the non square domain_in, such that
    % domain_in lies in the rectangle domain_out. Functions must be vectorized.
    %
    %
    % PUFun2DLS(f,domain_in,domain_out,'perf1',perf1,'pref2',pref2,..) constructs a
    % partition of unity approximation representing f, based on the options
    % passed into with varargin; that is PUFun(f,'perf1',perf1,'pref2',pref2,..)
    % is called. This preferences that can be set are:
    %
    % The max lengths of the patches before sampling is to occur:
    % 'MaxLengths', [d_1 d_2]
    %
    % *The degree indices from the standard degrees in each dimension for non
    % square domains : 'degreeIndex', [ind_1,ind_2].
    %
    % *The tolerance used for adaptation: 'tol', tol.
    %
    % *The degree indices from the standard degrees in each dimension for
    % square domains : 'ChebDegreeIndex', [ind_1,ind_2].
    %
    %
    % Here the degrees can be chosen from the set [3 5 9 17 33 65 129].
    % So if 'degreeIndex', [5 5 5], the max degree of any approximate will be
    % 33 in each direction.
    properties
        ChebRoot
        leafArray
        deg_in
        cheb_deg_in
        domain_in
        tol
        domain
        grid_opt = false;
    end
    
    methods
        
        function obj = PUFunLS(varargin)
            
            f = varargin{1};
            
            obj.domain_in = varargin{2};
            obj.domain = varargin{3};
            
            cheb_struct.domain_in = obj.domain_in;
            cheb_struct.domain = obj.domain;
            
            
            varargin(1:3) = [];
            args = varargin;
            
            while ( ~isempty(args) )
                if strcmpi(args{1}, 'degreeIndex')
                    cheb_struct.deg_in = args{2};
                elseif strcmpi(args{1}, 'ChebDegreeIndex')
                    cheb_struct.cheb_deg_in = args{2};
                elseif strcmpi(args{1}, 'MaxLengths')
                    cheb_struct.max_lengths = args{2};
                elseif strcmpi(args{1}, 'tol')
                    cheb_struct.tol = args{2};
                elseif strcmpi(args{1}, 'CourseDegreeIndex')
                    cheb_struct.cdeg_in = args{2};
                else
                    error(strcat(args{1},' is not a valid parameter.'));
                end
                args(1:2) = [];
            end
            
            dim = nargin(f);
            
            if dim==2
                obj.ChebRoot = LSPatch2D(cheb_struct);
            else
                obj.ChebRoot = LSPatch3D(cheb_struct);
            end
            
            obj.cheb_deg_in = obj.ChebRoot.cheb_deg_in;
            obj.deg_in = obj.ChebRoot.deg_in;
            obj.tol = obj.ChebRoot.tol;
            
            refine(obj,f);
            
            obj.ChebRoot.clean();
            
        end
        
        % diff_Tree = diff(obj,diff_dim,order)
        % This method computes thd approximation of the derivative
        %Input:
        %Output:
        %   diff_Tree  : PU approximation of derivative
        function diff_Tree = diff(obj,diff_dim,order)
            diff_Tree = copy(obj);
            
            for i=1:length(obj.leafArray)
                if isa(diff_Tree.leafArray{i},'LSPatch')
                    diff_Tree.leafArray{i}.coeffs = Diff(obj.leafArray{i},diff_dim,order);
                else
                    diff_Tree.leafArray{i}.values = evalfDiffGrid(obj.leafArray{i},diff_dim,order);
                end
            end
            
        end
        
        % refine(obj,f,grid_opt)
        % This method refines the tree by f(x).
        %Input:
        %   f          : the function to split on
        %   grid_opt   : boolean value indicating if
        %                function is evaluated for grids;
        %                must take cell array of grids
        function refine(obj,f,grid_opt)
            
            if nargin<3
                grid_opt = false;
            end
            
%             h = figure();
%             
%             plot(obj.domain_in); hold on; obj.ChebRoot.plotzone; hold off;
%             
%            
%             frame = getframe(h);
%             im = frame2im(frame);
%             [imind,cm] = rgb2ind(im,256);
%             
%             imwrite(imind,cm,'cool_mov2.gif','gif', 'Loopcount',inf); 
%             
%             close all
            
            while ~obj.ChebRoot.is_refined
                
                %split to get satisfactory covering before sampling
                %obj.Geomrefine();
                
                %then sample;
                Max = obj.ChebRoot.sample(f,grid_opt);
                
                 
                
                if obj.ChebRoot.is_leaf
                    obj.ChebRoot = obj.ChebRoot.splitleaf(Max);
                else
                    obj.ChebRoot.PUsplit(Max);
                end
                
            
%             h = figure();
%             
%             plot(obj.domain_in); hold on; obj.ChebRoot.plotzone; hold off;
%             
%             frame = getframe(h);
%             im = frame2im(frame);
%             [imind,cm] = rgb2ind(im,256);
%                 
%             imwrite(imind,cm,'cool_mov2.gif','gif','WriteMode','append'); 
%             
%             close all
            end
            
            
            if ~obj.ChebRoot.is_leaf
                obj.ChebRoot.findIndex([]);
                obj.leafArray = obj.ChebRoot.collectLeaves();
            else
                obj.leafArray = {obj.ChebRoot};
            end
            
            obj.ChebRoot.clean();
            
        end
        
                % refine(obj,f,grid_opt)
        % This method refines the tree by f(x).
        %Input:
        %   f          : the function to split on
        %   grid_opt   : boolean value indicating if
        %                function is evaluated for grids;
        %                must take cell array of grids
        function Geomrefine(obj)
            while ~obj.ChebRoot.is_geometric_refined
                if obj.ChebRoot.is_leaf
                    obj.ChebRoot = obj.ChebRoot.splitleafGeom();
                else
                    obj.ChebRoot.PUsplitGeom();
                end
            end
        end
    end
    
    
    methods(Access = protected)
        % Override copyElement method:
        function cpObj = copyElement(obj)
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            
            cpObj.ChebRoot = obj.ChebRoot.copy();
            
            cpObj.leafArray = cpObj.ChebRoot.collectLeaves();
        end
    end
    
end