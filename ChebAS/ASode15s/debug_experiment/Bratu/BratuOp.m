classdef BratuOp
    
    properties
        domain = [-1 1;-1 1];
        n
        disc
        leaves %references to leaves, so ok.
        lambda
        bound
        is_pack = false;
        vb
     end
    
    
    properties (Access=private)
        the_rho
        the_drho_dt
        the_period
        the_pC = 0.5
    end
    
    methods
        
        function r = BratuOp(leaves,lambda,is_pack)
            
            r.is_pack = is_pack;
            r.leaves = leaves;
            r.n = leaves{1}.degs;
            r.domain = leaves{1}.domain;
            
            r = setup_discretization(r);
            
            r.lambda = lambda;
            
            
            %r.vb = @(x,y,t) 3/4 + 1./(4*( 1+ exp(r.lambda*(-t-4*x+4*y)/32)));
        end
        
        
        
        function [U] = unpack(r,u)
            
            r.leaves{1}.Setvalues(u);
            U = reshape(r.leaves{1}.Getvalues(true),r.leaves{1}.degs);
            
            
        end
        
        function u = pack(r,U)
            
            % assume boundary information is within patch,
            % set at initialization
            
            ind = ~r.leaves{1}.outer_boundary;
            u = U(ind);
          
            
        end
        
        function u_t = timederiv(r,t,u)
             
            if r.is_pack
                r = r.setBoundary(t);
                [U] = r.unpack(u);
            else
                U = reshape(u,r.n);
            end
            
            Ux = r.disc.dx*U; Uxx = r.disc.dx*Ux; 
            Uy = U*r.disc.dy'; Uyy = Uy*r.disc.dy';
            
            u_t = (Uxx+Uyy) + r.lambda*exp(U);
            
            if r.is_pack
                u_t = r.pack(u_t);
            else
                boundary = r.leaves{1}.outer_boundary;
                u_t(boundary) = u(boundary);
                u_t = u_t(:);
            end
            
        end
        
        function J = jac(r,t,u)
             
            J = (r.disc.Dx2+r.disc.Dy2);
             
            ind = r.leaves{1}.outer_boundary;
            
            if r.is_pack
                J = J(~ind,~ind) + r.lambda*diag(exp(u));
            else
                J = J+r.lambda*diag(exp(u));
                I = eye(length(u));
                J(ind,:) = I(ind,:);
            end
             
        end
        
        function r = setBoundary(r,t)
            r.leaves{1}.SetBoundaryValues(@(x,y)r.bound(x,y,t));
        end
        
        function M = massMatrix(r)
            if r.is_pack
                M = eye(length(r.leaves{1}));
            else
           ind = r.leaves{1}.outer_boundary;
           Z = zeros(length(r.leaves{1}));
           M(ind,:) = Z(ind,:);
            end
        end
        
        function [U] = initial(r)
            
            r = setBoundary(r,0);
            
            r.leaves{1}.Setvalues(@(x,y)r.bound(x,y,0));
            U = r.leaves{1}.Getvalues();
            
        end
        
    end
    
    methods (Access=private)
        function r = setup_discretization(r)
            
            x = chebpts(r.n(1),r.domain(1,:));  
            y = chebpts(r.n(2),r.domain(2,:));
            
            % Discretization in [-1,1]^2
            dx = diffmat(r.n(1),1,r.domain(1,:));
            dy = diffmat(r.n(2),1,r.domain(2,:));
            
            Ix = eye(r.n(1)); Iy = eye(r.n(2));
            
            Dx = kron(Iy,dx);
            Dy = kron(dy,Ix);
            

            dx2 = dx^2;
            dy2 = dy^2;
            
            Dx2 = kron(Iy,dx2);
            Dy2 = kron(dy2,Ix);
            
            
            r.disc.x = x;
            r.disc.y = y;
            
            r.disc.dx = dx;
            r.disc.dy = dy;
            
            r.disc.dx2 = dx2;
            r.disc.dy2 = dy2;
            
            r.disc.Dx = Dx;
            r.disc.Dx2 = Dx2;
            
            r.disc.Dy = Dy;
            r.disc.Dy2 = Dy2;

        end
        
    end
    
    
    
    
end

