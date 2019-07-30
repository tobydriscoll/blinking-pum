classdef blinkone
    
    properties
        domain
        n 
		model
		map
		lid
		flux
		disc
	end

    
    properties (Dependent)
        pA, pS
		boundaryH, h_e
        rho
        drho_dt
		period
 		Qtop, Qleft
    end
    
    
    methods
        
        function r = blinkone(model,leaf,map,lid,flux)
			r.model = model;
			r.map = map;
			r.lid = lid;
			r.flux = flux;
            r.n = leaf.degs;
			r.domain = leaf.domain;           
            r = setup_discretization(r,leaf);
        end
        
        function rho = get.rho(r)
            rho = r.lid.rho;
        end
        function drho = get.drho_dt(r)
            drho = r.lid.drho_dt;
		end
		function T = get.period(r)
            T = r.lid.period;
		end
		function A = get.pA(r)
            A = r.model.A;
        end
		function S = get.pS(r)
            S = r.model.S;
		end 
		function x = get.h_e(r)
            x = r.model.h_slideover;
		end
 		function x = get.boundaryH(r)
            x = r.model.h_boundary;
		end
		function x = get.Qtop(r)
            x = r.flux.Qtop;
		end
		function x = get.Qleft(r)
            x = r.flux.Qleft;
		end
 
        function massmat = massMatrix(r)
            MH = ones(r.disc.num.h,1);   MP = zeros(r.disc.num.p,1);
            N = r.disc.num.h + r.disc.num.p;
            massmat = spdiags([MH;MP],0,N,N);
		end                                             
                       
    end
    
    methods (Hidden=true)
        
        %% Functions for computing the solution.
        function [W,Ws,J,Js] = grad(r,t,V,factor)
            % Gradient of a scalar function (physical domain and strip)
            if nargin < 4
                [~,~,factor] = r.disc.strip.grid(t);
            end
            D = r.disc.strip.dm;
            Ws = D.x*V + 1i*V*D.y(t)';
            W = conj(factor).*Ws;
            if nargout > 2
                % Jacobian
                I = r.disc.eye;
                Js = kron(I.y,D.x) + 1i*kron(D.y(t),I.x);
                J = conj(r.disc.Diag(factor))*Js;
            end
        end
        
        function W = div(r,t,V,factor)
            % Divergence of a vector field (physical domain)
            if nargin < 4
                [~,~,factor] = r.disc.strip.grid(t);
            end
            D = r.disc.strip.dm;
            W = real( factor.*(D.x*V - 1i*V*D.y(t)') );
        end
        
        function [H,P] = unpack(r,u)
            % Break ODE unknown into h and p function values
            H = repmat(r.boundaryH,r.n);
            H(~r.disc.boundary.loc_outer) = u(1:r.disc.num.h);
            P = r.disc.unvec(u(r.disc.num.h+1:end));
        end
        
        function u = pack(r,H,P)
            u = [r.disc.vec(H(~r.disc.boundary.loc_outer));r.disc.vec(P)];
        end
        
        function u_t = timederiv(r,t,u)
            % Compute f in M*u_t = f(t,u)
            [H,P] = unpack(r,u);
            [XS,YS,factor] = r.disc.strip.grid(t);
            xc = r.disc.square.points.x;  yc = r.disc.square.points.y;
            [gradH,gradHs] = r.grad(t,H,factor);
            [gradP,gradPs] = r.grad(t,P,factor);
            
            Psi = H.^3/3;
            Q = -Psi.*gradP;
            DivQ = r.div(t,Q,factor);
            DivH = r.div(t,gradH,factor);
            
            Motion = (r.disc.strip.dydt(t)'/r.map.strip.dydyinv(t)) .* imag(gradHs);
            H_t = -DivQ - Motion;
            P_t = P + r.pS*DivH + r.pA*H.^-3;
            
            % Flux condition on left side (drainage)
            xleft = XS(1,:);  yleft = YS(1,:);
            absfp = r.map.eye.absdzdz(xleft,yleft);
            Qdotn = -Psi(1,:).*real(gradPs(1,:))./absfp;
            Qout = r.Qleft(t,yc);
            P_t(1,:) = Qdotn - Qout(:)';  % east
            
            % No flux on two sides
            xb = XS(end,:);  yb = YS(end,:);
            absfp = r.map.eye.absdzdz(xb,yb);
            Qdotn = Psi(end,:).*real(gradPs(end,:))./absfp;
            P_t(r.n(1),:) = Qdotn;   % west

            xb = XS(:,1);  yb = YS(:,1);
            absfp = r.map.eye.absdzdz(xb,yb);
            Qdotn = -Psi(:,1).*imag(gradPs(:,1))./absfp;
            P_t(:,1) = Qdotn;                % south
            
            % Flux on moving side
            xtop = XS(:,end);  ytop = YS(:,end);
            absfp = r.map.eye.absdzdz(xtop,ytop);
            Qdotn = Psi(:,r.n(2)).*imag(gradPs(:,r.n(2)))./absfp;
            Qin = r.Qtop(t,xc);
            P_t(:,r.n(2)) = Qdotn - Qin + r.drho_dt(t).*absfp.*(H(:,r.n(2))-r.h_e/2);
            
            u_t = r.pack(H_t,P_t); 
           % fprintf('\b\b\b\b\b\b%.4f',t)

        end
        
        function J = jac(r,t,u)
            [H,P] = r.unpack(u);
            [XS,YS,factor] = r.disc.strip.grid(t);
%            Factor = r.disc.Diag(factor);
            
            [~,~,JH,Jgrads] = r.grad(t,H,factor);
            [gradP,~,JP] = r.grad(t,P,factor);
            
            Psi = H.^3/3;
            dPsi_dh = H.^2;
            
            %Q = -Psi.*gradP;
            DiagPsi = r.disc.Diag(Psi);
            
            i1r = 1:numel(H);
            i1c = 1:numel(H);
            i2r = i1r(end) + (1:numel(H));
            i2c = i1c(end) + (1:numel(H));
            J = zeros(i2r(end),i2c(end));
            
            %DivQ = div(t,Q,factor);
            cJs = conj(Jgrads); 
            B = full(factor(:).*cJs);
            Br = real(B);  Bi = imag(B);
            A = dPsi_dh(:).*gradP(:);
            J(i1r,i1c) = -Br.*real(A)' + Bi.*imag(A)';
            A = Psi(:).*full(JP);
            J(i1r,i2c) = Br*real(A) - Bi*imag(A);

            %Motion = bsxfun(@times, dyc_dt(t)'./dyc_dys(t), imag(gradHs) );
            v = r.disc.strip.dydt(t)/r.map.strip.dydyinv(t);
            JMo = kron( spdiags(v,0,r.n(2),r.n(2))*r.disc.strip.dm.y(t), r.disc.eye.x );
            %H_t = -DivQ - Motion;
            J(i1r,i1c) = J(i1r,i1c) + JMo;
            
            %DivH = div(t,gradH,factor);
            %P_t = P + pS*DivH + pA*H.^-3;
            A = full(JH);
            J(i2r,i1c) = r.pS*(Br*real(A)-Bi*imag(A)) - 3*r.pA*r.disc.Diag(H.^(-4));
            J(i2r,i2c) = r.disc.eye.all;
            
            % Flux conditions
            %P_t([1 n(1)],:) = real(gradPs([1 n(1)],:));  % east and west
            %P_t(:,1) = imag(gradPs(:,1));                % south
            Jxs = real(Jgrads);  Jys = imag(Jgrads);
            bdy = r.disc.boundary;
            
            Jn = DiagPsi*Jxs;
            xleft = XS(1,:);  yleft = YS(1,:);
            absfp = r.map.eye.absdzdz(xleft,yleft);
            J(i2r(bdy.E),:) = [0*Jxs(bdy.E,:) -diag(1./absfp)*Jn(bdy.E,:)];
            
            xb = XS(end,:);  yb = YS(end,:);
            absfp = r.map.eye.absdzdz(xb,yb);
            J(i2r(bdy.W),:) = [0*Jxs(bdy.W,:) diag(1./absfp)*Jn(bdy.W,:)];
            
            Jn = DiagPsi*Jys;
            xb = XS(:,1);  yb = YS(:,1);
            absfp = r.map.eye.absdzdz(xb,yb);
            J(i2r(bdy.S),:) = [0*Jys(bdy.S,:) -diag(1./absfp)*Jn(bdy.S,:)];
            
            % Moving side
            xtop = XS(:,end);  ytop = YS(:,end);
            absfp = r.map.eye.absdzdz(xtop,ytop);
            %Qdotn = Psi(:,n(2)).*imag(gradPs(:,n(2)))./absfp;
            %P_t(:,n(2)) = Qdotn + drho_dt(t).*absfp.*(H(:,n(2))-He);
            J(i2r(bdy.N),:) = [diag(absfp)*r.drho_dt(t)*r.disc.eye.all(bdy.N,:),...
                diag(1./absfp)*Jn(bdy.N,:) ];
            %Qdotn = Pphi(:,n(2)).*Phi_ys(:,n(2))./amp;
            %Phi_t(:,n(2)) = Qdotn + drho_dt(t).*amp.*Phi(:,n(2));
            
            %J=full([J1;J2]);
            
            % Boundary h values do not appear
            idx = find(r.disc.boundary.loc_outer);
            J(idx,:) = [ ];    J(:,idx) = [ ];
            
        end
              
    end
    
    methods (Access=private)
        function r = setup_discretization(r,leaf)
            
            % Discretization in [-1,1]^2
            xc = chebpts(r.n(1),r.domain(1,:));  Dxc = diffmat(r.n(1),1,r.domain(1,:));
            yc = chebpts(r.n(2),r.domain(2,:));  Dyc = diffmat(r.n(2),1,r.domain(2,:));
            [XC,YC] = ndgrid(xc,yc);
            Ix = speye(r.n(1));  Iy = speye(r.n(2));
            I = speye(prod(r.n));
			
			% Chain rule terms
			map = r.map;            
            dxs_dxc = map.strip.dxdx(xc);
            Dxs = diag(1./dxs_dxc)*Dxc;
            Dys = @(t) Dyc*map.strip.dydyinv(t);  % leaves out the d/dt term for below
            
            function [XS,YS,factor] = stripgrid(t,xc_o,yc_o)               
                if nargin==1
					[XS,YS] = ndgrid(map.strip.x(xc),map.strip.y(t,yc));
				else
					[XS,YS] = ndgrid(map.strip.x(xc_o),map.strip.y(t,yc_o));
				end
				if nargout > 2
					factor = 1./map.eye.dzdz(XS,YS);
				end
			end
                      
            % Identify the boundaries in the grid
            east = false(r.n);   east(1,:) = true;
            west = false(r.n);   west(r.n(1),:) = true;
            south = false(r.n);  south(:,1) = true;
            north = false(r.n);  north(:,r.n(2)) = true;
            
			% store it all
            discr.eye = struct('x',Ix,'y',Iy,'all',I);
            discr.vec = @(U) U(:);  discr.unvec = @(u) reshape(u,r.n);
            discr.Diag = @(U) spdiags(discr.vec(U),0,prod(r.n),prod(r.n));
            discr.num = struct('h',sum(~leaf.outer_boundary),'p',prod(r.n));
            discr.square.points = struct('x',xc,'y',yc);
            discr.square.grid = struct('X',XC,'Y',YC);
            discr.square.dm = struct('x',Dxc,'y',Dyc);
            discr.strip.points = struct('x',map.strip.x(xc),'y',@(t) map.strip.y(t,yc));
            discr.strip.grid = @stripgrid;
            discr.strip.dydt = @(t) map.strip.dydt(t,yc);
            discr.strip.dm = struct('x',Dxs,'y',Dys);
            discr.boundary = struct('N',north,'S',south,'E',east,'W',west);
            discr.boundary.all = east | west | north | south;
            discr.boundary.loc_outer = leaf.outer_boundary;           
            r.disc = discr;
        end
        
        function v = volume_vals(r,t,U,wx,wy)
            if nargin < 5
                [~,wx] = chebpts(r.n(1));  [~,wy] = chebpts(r.n(2));
            end
            % Jacobian for comp. to strip map
            A = (r.disc.strip.map.deriv.x*(r.rho(t)+1)/2) .* U;
            [XS,YS] = r.disc.strip.grid(t);
            % Jacobian for strip to eye map
            A = A.*r.disc.eye.map.absderiv(XS,YS).^2;
            v = (wx*A*wy');
        end
                        
    end
    
    
    
    
end