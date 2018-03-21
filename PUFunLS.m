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
    % *The degree indices from the standard degrees in each dimension for
    % square domains : 'ChebDegreeIndex', [ind_1,ind_2].
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
        domain_in
        tol
        domain
        grid_opt = false;
        movie;
    end
    
    methods
        
        function obj = PUFunLS(varargin)
            
            if length(varargin)==1
                varargin = varargin{:};
            end
            
            if isstruct(varargin)
                
                obj.grid_opt = varargin.grid_opt;
                f = varargin.op;
                obj.domain = varagin.domain;
                obj.domain_in = varagin.domain_in;
                obj.deg_in = varagin.deg_in;
                obj.cheb_deg_in = varagin.cheb_deg_in;
                obj.tol = varagin.tol;
                
                dim = nargin(f);
                
                if dim==2
                    obj.ChebRoot = LSPatch2D('domain',obj.domain_in,'boundingbox',obj.domain,'degreeIndex',obj.deg_in,'ChebDegreeIndex',obj.cheb_deg_in,'tol',obj.tol);
                else
                    obj.ChebRoot = LSPatch3D('domain',obj.domain_in,'boundingbox',obj.domain,'degreeIndex',obj.deg_in,'ChebDegreeIndex',obj.cheb_deg_in,'tol',obj.tol);
                end
                
            else
                if length(varargin)==3
                    f = varargin{1};
                    
                    obj.domain_in = varargin{2};
                    obj.domain = varargin{3};
                    
                    dim = nargin(f);
                    
                    if dim==2
                        obj.ChebRoot = LSPatch2D('InnerDomain',obj.domain_in,'domain',obj.domain);
                    else
                        obj.ChebRoot = LSPatch3D('InnerDomain',obj.domain_in,'domain',obj.domain);
                    end
                else
                    f = varargin{1};
                    
                    obj.domain_in = varargin{2};
                    obj.domain = varargin{3};
                    varargin(1:3) = [];
                    args = varargin;
                    while ( ~isempty(args) )
                        if strcmpi(args{1}, 'degreeIndex')
                            obj.deg_in = args{2};
                        elseif strcmpi(args{1}, 'ChebDegreeIndex')
                            obj.cheb_deg_in = args{2};
                        end
                        args(1:2) = [];
                    end
                    
                    
                    varargin = {varargin{:},'InnerDomain',obj.domain_in,'domain',obj.domain};
                    
                    dim = nargin(f);
                    
                    if dim==2
                        obj.ChebRoot = LSPatch2D(varargin);
                    else
                        obj.ChebRoot = LSPatch3D(varargin);
                    end
                    
                end
            end
            
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
            
            h = figure();
            
            obj.ChebRoot.plotzone; hold on; plot(obj.domain_in);
            
           
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            
            imwrite(imind,cm,'cool_mov.gif','gif', 'Loopcount',inf); 
            
            close all
            
            while ~obj.ChebRoot.is_refined
                
                %split to get satisfactory covering before sampling
                %obj.Geomrefine();
                
                %then sample;
                Max = obj.ChebRoot.sample(f,grid_opt);
                
                obj.ChebRoot.plotzone; hold on; plot(obj.domain_in);
                
                frame = getframe(gca); 
                im = frame2im(frame); 
                [imind,cm] = rgb2ind(im,256); 
                
                if obj.ChebRoot.is_leaf
                    obj.ChebRoot = obj.ChebRoot.splitleaf(Max);
                else
                    obj.ChebRoot.PUsplit(Max);
                end
                
            
            h = figure();
            
            obj.ChebRoot.plotzone; hold on; plot(obj.domain_in); hold off;
            
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
                
            imwrite(imind,cm,'cool_mov.gif','gif','WriteMode','append'); 
            
            close all
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