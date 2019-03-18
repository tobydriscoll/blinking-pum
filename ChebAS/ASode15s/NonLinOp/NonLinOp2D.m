classdef NonLinOp2D < NonLinOp
    %NonLinOp Abstract class used for our PDE solvers on rectangular
    %domains.
    
    properties
        
        
        Chebs %cell array of ChebPatches
        
        boundf = [];
        
        rfun
        
        jfun
        
    end
    
    properties(Access = protected)
        
        disc_f %full discritization
        
        disc_c %coarse discritization
    end
    
    methods
        function obj = NonLinOp2D(Chebs,rfun,jfun,boundf)
            
            if ~iscell(Chebs)
                obj.Chebs = {Chebs};
            else
                obj.Chebs = Chebs;
            end
            
            if nargin==2
                if ~iscell(boundf)
                    obj.boundf = {boundf};
                else
                    obj.boundf = {boundf};
                end
                
            end
            
            obj.disc_f = SetUpDiscritization(obj);
            
            for i=1:length(Chebs)
                Chebs{i}.Coarsen();
            end
            
            obj.disc_c = SetUpDiscritization(obj);
            
            for i=1:length(Chebs)
                Chebs{i}.Refine();
            end
            
            obj.rfun = rfun;
            
            obj.jfun = jfun;
            
            obj.boundf = {boundf};
        end
        
        function R = residual(obj,u,params)
           
            
            SetBoundaryVals(obj);
            
            u = unpack(obj,u);
            
            R = obj.rfun(obj,u);
            
            R = pack(obj,R);
           
        end
        
        function J = jac(obj,u,params)
            
            SetBoundaryVals(obj)
            
            u = unpack(obj,u);
            
            J = obj.jfun(obj,u);
            
            ind = true(0,1);
            
            for i=1:length(obj.Chebs)
                if obj.Chebs{i}.is_packed
                    ind = [ind;~obj.Chebs{i}.outer_boundary];
                else
                    ind = [ind;true(length(obj.Chebs{i}),1)];
                end
            end
            
            J = J(ind,ind);
            
        end
        
        function [U,Ux,Uxx,Uy,Uyy] = computeDerives(obj,u)
            
            D = GetDisc(obj);
            
            U = reshape(u,obj.Chebs{1}.degs);
            
            Ux = D.dx*U; Uxx = D.dx*Ux;
            Uy = U*D.dy'; Uyy = Uy*D.dy';
            
        end
        
        function D = GetDisc(obj)
            
            if obj.Chebs{1}.iscoarse
                D = obj.disc_c;
            else
                D = obj.disc_f;
            end
            
        end
        
        
        
        function [v] = pack(obj,u)
            
            v = cell(length(obj.Chebs),1);
            
            index = 0;
            
            for i=1:length(obj.Chebs)
                
                sub_ind = index + (1:prod(obj.Chebs{i}.degs));
                
                u_i = u(sub_ind);
                
                if obj.Chebs{i}.is_packed
                    v{i} = u_i(~obj.Chebs{i}.outer_boundary);
                else
                    v{i} = u_i;
                end
                
                index = index + prod(obj.Chebs{i}.degs);
                
            end
            
            v = cell2mat(v);
            
        end
        
        function [v] = unpack(obj,u)
            
            v = cell(length(obj.Chebs),1);
            
            index = 0;
            
            for i=1:length(obj.Chebs)
                
                sub_ind = index + (1:length(obj.Chebs{i}));
                
                u_i = u(sub_ind);
                
                if obj.Chebs{i}.is_packed
                    obj.Chebs{i}.Setvalues(u_i);
                    v{i} = obj.Chebs{i}.Getvalues(true);
                else
                    v{i} = u_i;
                end
                
                index = index + length(obj.Chebs{i});
                
            end
            
            v = cell2mat(v);
            
        end
        
        function SetBoundaryVals(obj)
            
            for i=1:length(obj.Chebs)
                if obj.Chebs{i}.is_packed
                    obj.Chebs{i}.SetBoundaryValues(obj.boundf{i});
                end
            end
            
        end
        
        
    end
    
    methods (Access  = protected)
        function D = SetUpDiscritization(obj)
            
            degs = obj.Chebs{1}.degs();
            
            domain = obj.Chebs{1}.domain;
            
            
            D.x = chebpts(degs(1),domain(1,:));
            D.y = chebpts(degs(2),domain(2,:));
            
            D.dx = diffmat(degs(1),1,domain(1,:));
            D.dy = diffmat(degs(2),1,domain(2,:));
            
            D.Ix = eye(degs(1)); D.Iy = eye(degs(2));
            
            D.Dx = kron(D.Iy,D.dx);
            D.Dy = kron(D.dy,D.Ix);
            
            
            D.dx2 = D.dx^2;
            D.dy2 = D.dy^2;
            
            D.Dx2 = kron(D.Iy,D.dx2);
            D.Dy2 = kron(D.dy2,D.Ix);
        end
    end 
        
    
end



