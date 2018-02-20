classdef PUFun2DLS < handle
    
    properties
        ChebRoot
        TreeGrid
        leafArray
        Errs
        Nums
        deg_in
        tol
        domain
    end
    
    methods
        
        function obj = PUFun2DLS(domain,in_domain,deg_in,max_lengths,f,tol,grid_opt,ChebRoot)
            obj.domain = domain;
            obj.deg_in = deg_in;
            [dim,~] = size(domain);
            obj.tol = tol;
            
            if nargin < 7
                grid_opt = false;
                obj.ChebRoot = LSPatch2D(in_domain,max_lengths,domain,domain,domain);
            elseif nargin < 8
                obj.ChebRoot = LSPatch2D(in_domain,max_lengths,domain,domain,domain);
            else
                obj.ChebRoot = ChebRoot;
            end
            
            if nargin<7
                refine(obj,f,grid_opt);
            else
                obj.TreeGrid = obj.ChebRoot.leafGrids();
                
                obj.leafArray = obj.ChebRoot.collectLeaves();
                
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
            
            while ~obj.ChebRoot.is_refined
                
                obj.ChebRoot.sample(f,grid_opt);
                
                if obj.ChebRoot.is_leaf
                    obj.ChebRoot = obj.ChebRoot.splitleaf();
                else
                    obj.ChebRoot.PUsplit(inf);
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