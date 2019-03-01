classdef BurgersOp
    
    properties
        domain = [-1 1;-1 1];
        n
        disc
        leaves %references to leaves, so ok.
        R
        ub
        vb
        is_pack
        
     end
    
    
    properties (Access=private)
        the_rho
        the_drho_dt
        the_period
        the_pC = 0.5
    end
    
    methods
        
        function r = BurgersOp(leaves,R,is_pack)

            r.is_pack = is_pack;
            r.leaves = leaves;
            r.n = leaves{1}.degs;
            r.domain = leaves{1}.domain;
            
            r = setup_discretization(r);
            
            r.R = R;
            
            r.ub = @(x,y,t) 3/4 - 1./(4*( 1+ exp(r.R*(-t-4*x+4*y)/32)));
            r.vb = @(x,y,t) 3/4 + 1./(4*( 1+ exp(r.R*(-t-4*x+4*y)/32)));
            
            %r.ub = @(x,y,t) ones(size(x));
            %r.vb = @(x,y,t) ones(size(x));
        end
        
        
        
        function [U,V] = unpack(r,u)
            
            len = length(r.leaves{1});
            
            r.leaves{1}.Setvalues(u(1:len));
            U = reshape(r.leaves{1}.Getvalues(true),r.leaves{1}.degs);
            
            r.leaves{2}.Setvalues(u(len+1:end));
            V = reshape(r.leaves{2}.Getvalues(true),r.leaves{2}.degs);
            
        end
        
        function u = pack(r,U,V)
            
            % assume boundary information is within patch,
            % set at initialization
            
            ind = ~r.leaves{1}.outer_boundary;
            u = [U(ind);V(ind)];
          
            
        end
        
        function u_t = timederiv(r,t,y)
             
            if r.is_pack
             r = r.setBoundary(t);
             
             [U,V] = r.unpack(y);
             
            else
             degs = r.leaves{1}.degs;
             U = zeros(degs);
             V = zeros(degs);
             Len = prod(degs);
             U(:) = y(1:Len); %x_velocity
             V(:) = y(Len+(1:Len)); %y_velocity
            end
             
             Ux = r.disc.dx*U; Uxx = r.disc.dx*Ux; 
             Uy = U*r.disc.dy'; Uyy = Uy*r.disc.dy';
             
             Vx = r.disc.dx*V; Vxx = r.disc.dx*Vx; 
             Vy = V*r.disc.dy'; Vyy = Vy*r.disc.dy';
             
             F1 = (1/r.R)*(Uxx+Uyy) - (U.*Ux+V.*Uy);
             F2 = (1/r.R)*(Vxx+Vyy) - (U.*Vx+V.*Vy);
             
             border = r.leaves{1}.outer_boundary;
    
             P = r.leaves{1}.points();
             
             F1(border) = U(border) - r.ub(P(border,1),P(border,2),t);
             F2(border) = V(border) - r.vb(P(border,1),P(border,2),t);
             
             if r.is_pack
                u_t = r.pack(F1,F2);
             else
                u_t = [F1(:);F2(:)];
             end

             
        end
        
        function J = jac(r,t,y)
            
            if r.is_pack
            r = r.setBoundary(t);
            
            [U,V] = r.unpack(y);
            else
            
             degs = r.leaves{1}.degs;
             U = zeros(degs);
             V = zeros(degs);
             Len = prod(degs);
             U(:) = y(1:Len); %x_velocity
             V(:) = y(Len+(1:Len)); %y_velocity
            
            end
            Ux = r.disc.dx*U; Uy = U*r.disc.dy';
            Vx = r.disc.dx*V; Vy = V*r.disc.dy';
            
            U = U(:); Ux = Ux(:); Uy = Uy(:);
            V = V(:); Vx = Vx(:); Vy = Vy(:);
            
            border = r.leaves{1}.outer_boundary;
            
            I = eye(prod(r.leaves{1}.degs));
            Z = zeros(prod(r.leaves{1}.degs));
            
            J11 = (1/r.R)*(r.disc.Dx2+r.disc.Dy2) - diag(U)*r.disc.Dx - diag(Ux)-diag(V)*r.disc.Dy;
            J12 = -diag(Uy);
            J11(border,:) = I(border,:); J12(border,:) = Z(border,:);
            
            J21 = -diag(Vx);
            J22 = (1/r.R)*(r.disc.Dx2+r.disc.Dy2) - diag(U)*r.disc.Dx - diag(V)*r.disc.Dy - diag(Vy);
            J21(border,:) = Z(border,:); J22(border,:) = I(border,:);
            
            if r.is_pack
                ind1 = ~r.leaves{2}.outer_boundary;
                
                J11 = J11(ind1,ind1);
                J12 = J12(ind1,ind1);

                ind2 = ~r.leaves{2}.outer_boundary;

                J21 = J21(ind2,ind2);
                J22 = J22(ind2,ind2);

            end
            
           J = [J11 J12;J21 J22];
        end
        
        function r = setBoundary(r,t)
            r.leaves{1}.SetBoundaryValues(@(x,y)r.ub(x,y,t));
            r.leaves{2}.SetBoundaryValues(@(x,y)r.vb(x,y,t));
        end
        
        function M = massMatrix(r)
            if r.is_pack
                M = eye(length(r.leaves{1})+length(r.leaves{2}));
            else
                I = eye(prod(r.leaves{1}.degs));
                border = r.leaves{1}.outer_boundary;
                I(border,:) = 0;
                M = blkdiag(I,I);
            end

        end
        
        function [U,V] = initial(r)
            
            r = setBoundary(r,0);
            
            r.leaves{1}.Setvalues(@(x,y)r.ub(x,y,0));
            U = r.leaves{1}.Getvalues();
            
            r.leaves{2}.Setvalues(@(x,y)r.vb(x,y,0));
            V = r.leaves{2}.Getvalues();
        end
        
    end
    
    methods (Access=private)
        function r = setup_discretization(r)
            
            degs = r.n;
            
            domain = r.domain;
            
            
            discr.x = chebpts(degs(1),domain(1,:));
            discr.y = chebpts(degs(2),domain(2,:));
            
            discr.dx = diffmat(degs(1),1,domain(1,:));
            discr.dy = diffmat(degs(2),1,domain(2,:));
            
            discr.Ix = eye(degs(1)); discr.Iy = eye(degs(2));
            
            discr.Dx = kron(discr.Iy,discr.dx);
            discr.Dy = kron(discr.dy,discr.Ix);
            
            
            discr.dx2 = discr.dx^2;
            discr.dy2 = discr.dy^2;
            
            discr.Dx2 = kron(discr.Iy,discr.dx2);
            discr.Dy2 = kron(discr.dy2,discr.Ix);
            
            r.disc = discr;

        end
        
    end
    
    
    
    
end