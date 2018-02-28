classdef PUFun2DLS < handle
    
    properties
        ChebRoot
        TreeGrid
        leafArray
        Errs
        Nums
        deg_in
        cheb_deg_in
        tol
        domain
        grid_opt = false;
    end
    
    methods
        
        %function obj = PUFun2DLS(domain,in_domain,deg_in,cheb_deg_in,max_lengths,f,tol,grid_opt,ChebRoot)
        function obj = PUFun2DLS(varargin)
            
            if length(varargin)==1
                varargin = varargin{:};
            end
            
            if isstruct(varargin)
                obj.grid_opt = varargin.grid_opt;
                f = varargin.op;
                obj.domain = varagin.domain;
                obj.deg_in = varagin.deg_in;
                obj.tol = varagin.tol;
                
                obj.ChebRoot = ChebPatch('domain',obj.domain,'degreeIndex',obj.deg_in,'tol',obj.tol);
                
            else
                if length(varargin)==1
                    f = varargin;
                    obj.domain = repmat([-1 1],nargin(f),1);
                    obj.ChebRoot = ChebPatch('domain',obj.domain);
                    
                elseif length(varargin)==2
                    f = varargin{1};
                    obj.domain = varargin{2};
                    obj.ChebRoot = ChebPatch('domain',obj.domain);
                    
                else
                    f = varargin{1};
                    varargin(1) = [];
                    args = varargin;
                    while ( ~isempty(args) )
                        if strcmpi(args{1}, 'gridOption')
                            obj.grid_opt = args{2};
                        elseif strcmpi(args{1}, 'domain')
                            obj.domain = args{2};
                        elseif strcmpi(args{1}, 'degreeIndex')
                            obj.deg_in = args{2};
                        end
                        args(1:2) = [];
                    end
                    
                    
                    if isempty(obj.domain)
                        obj.domain = repmat([-1 1],nargin(f),1);
                        varargin = {varargin{:} , 'domain',obj.domain};
                    end
                    
                    obj.ChebRoot = ChebPatch(varargin);
                    
                end
            end
            
            refine(obj,f);
            
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
            
            
            while ~obj.ChebRoot.is_refined
                
                obj.ChebRoot.IsGeometricallyRefined();
                Max = obj.ChebRoot.sample(f,grid_opt);
                
                if obj.ChebRoot.is_leaf
                    obj.ChebRoot = obj.ChebRoot.splitleaf(Max);
                else
                    obj.ChebRoot.PUsplit(Max);
                end
            end
            
            
            if ~obj.ChebRoot.is_leaf
                obj.ChebRoot.findIndex([]);
                obj.leafArray = obj.ChebRoot.collectLeaves();
            else
                obj.leafArray = {obj.ChebRoot};
            end
            
            obj.ChebRoot.clean();
            
        end
    end
end