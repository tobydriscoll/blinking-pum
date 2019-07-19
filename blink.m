classdef blink
    
    properties
        n = [19 25]
        tspan = [0 5.258]
        pA = 0
        pS = 1e-3
        boundaryH = 13
        h_e = 0
        
        initcond = 'laplace'
        initvolume = 12;
        
        drainvolume = 0;
        supplyvolume = 0;
        Qtop
        Qleft
        
        disc
        
        odetol = 1e-4
        solvetime
        solution

        domain = [-1 1;-1 1];

        finalstate
    end

    
    properties (Dependent)
        p, h
        rho
        drho_dt
        period
        percentClosed
        Hmax, Hmin, Pmax, Pmin
    end
    
    properties (Access=private)
        the_rho
        the_drho_dt
        the_period
        the_pC = 0.5
    end
    
    methods
        
        function r = blink(pctClosed,degs,domain,boundaryH,drainvolume)
            % Require Chebfun
            assert(exist('chebpts','file')>0,'Add Chebfun to path');
            
            if nargin > 0, r.percentClosed = pctClosed;  end
            
            if nargin > 1, r.n = degs; end
            
            if nargin > 2, r.domain = domain; end
            
            if nargin > 3, r.boundaryH = boundaryH; end
            
            if nargin > 4, r.drainvolume = drainvolume; r.supplyvolume = drainvolume; end
            
            if nargin>0
                r = setup_discretization(r);
            
                [r.Qtop,r.Qleft] = fluxfuns(r);
            end
            
        end
        
        function r = set.percentClosed(r,pc)
            [r.the_rho,r.the_drho_dt,r.the_period] = lidmotion_realistic(pc);
            r.the_pC = pc;
        end
        function pc = get.percentClosed(r)
            pc = r.the_pC;
        end
        function rho = get.rho(r)
            rho = r.the_rho;
        end
        function drho = get.drho_dt(r)
            drho = r.the_drho_dt;
        end
        function T = get.period(r)
            T = r.the_period;
        end
        
        function massmat = massMatrix(r)
            MH = ones(r.disc.num.h,1);   MP = zeros(r.disc.num.p,1);
            N = r.disc.num.h + r.disc.num.p;
            massmat = spdiags([MH;MP],0,N,N);
        end
        
        function r = solve(r,tspan)
            r = setup_discretization(r);
            r.tspan = tspan;
            [r.Qtop,r.Qleft] = fluxfuns(r);
            
            %% options
            MH = ones(r.n-2);   MP = zeros(r.n);
            N = r.disc.num.h + r.disc.num.p;
            massmat = spdiags([MH(:);MP(:)],0,N,N);
            opt = odeset('mass',massmat,'reltol',r.odetol,'abstol',r.odetol);
            opt = odeset(opt,'jacobian',@(t,y) jac(r,t,y));
            opt = odeset(opt,'outputfcn',@(t,y,f) ode_status(r,t,y,f));

            %% Initial condition
            [H,P,flag] = r.initial;
            if isempty(flag)
                u0 = r.pack(H,P);
            else
                u0 = H; 
                opt = odeset(opt,'initialslope',P);
            end
            
            %% solve
            fprintf('\n\nt = 0.0000')
            tic
            sol = ode15s(@(t,y) timederiv(r,t,y),r.tspan,u0,opt);
            r.solvetime = toc;
            fprintf('\n\n')
            
            %% extract final state (for solution continuation)
            [ufinal,dufinal] = deval(sol,sol.x(end));
            [H,P] = unpack(r,ufinal);
            r.finalstate.H = chebfun2(H'); r.finalstate.P = chebfun2(P');
            [dH,dP] = unpack(r,dufinal); 
            dH([1 end],:) = 0;  dH(:,[1 end]) = 0;
            r.finalstate.dH = chebfun2(dH'); r.finalstate.dP = chebfun2(dP');
            
            %% save memory
            if r.odetol > 1e-7
                sol.y = single(sol.y);
                sol.idata.dif3d = single(sol.idata.dif3d);
            end
            r.solution = sol;
            
            % may not have finished, and rounding to single can "change"
            % the final time
            r.tspan = [sol.x(1) sol.x(end)];
            
        end
        
        function v = get.Hmin(r)
            % Minimum value over the collocation values (not the
            % interpolant).
            v = min(r.boundaryH,min(min( r.solution.y(1:r.disc.num.h,:) )));
        end
        
        function v = get.Hmax(r)
            % Maximum value over the collocation values (not the
            % interpolant).
            v = max(r.boundaryH,max(max( r.solution.y(1:r.disc.num.h,:) )));
        end
        
        function v = get.Pmin(r)
            % Minimum value over the collocation values (not the
            % interpolant).
            v = min(min( r.solution.y(1+r.disc.num.h:end,:) ));
        end
        
        function v = get.Pmax(r)
            % Maximum value over the collocation values (not the
            % interpolant).
            v = max(max( r.solution.y(1+r.disc.num.h:end,:) ));
        end
        
        
        %% Evalauting the numerical solution.
        
        % Get a nice vector of times (equally spaced within the three
        % phases). This is hard-coded to go with the lidmotion function.
        function t = times(r,m)
            if nargin < 2, m = 240; end
            if length(m)==3
                m1 = m(1); m2 = m(2); m3 = m(3);
            else
                m1 = ceil(0.35*m);
                m3 = ceil(0.3*m);
                m2 = m - m1 - m3;
            end
            
            % times in one blink cycle
            t1 = linspace(0,0.19,m1+1);
            t2 = linspace(0.19,5.176,m2+1);
            t3 = linspace(5.176,5.258,m3);
            t = [ t1(1:m1) t2(1:m2) t3 ]';
            
            % extend to multiple cycles
            m = length(t);
            while t(end) < r.tspan(2)-1e-8
                t = [t; t(end-m+2:end)+r.period];
            end
            
            % clip to the solution interval
            t(t>r.tspan(2)) = [];
        end
        
        function [X,Y] = grid(r,t)
            g = r.disc.square.grid;
            [X,Y] = r.disc.map(t,g.X,g.Y);
        end
        
        function [x,y] = shape(r,t)
            [X,Y] = r.grid(t);
            x = chebfun2(X);
            y = chebfun2(Y);
        end
        
        function H = evalH(r,t)
            H = repmat(r.boundaryH,r.n);
            H(~r.disc.boundary.all) = deval(r.solution,t,1:r.disc.num.h);
        end
        
        function P = evalP(r,t)
            num = r.disc.num;
            P = r.disc.unvec(deval(r.solution,t,num.h+(1:num.p)));
        end
        
        function f = get.h(r)
            f = @eval;
            function h = eval(t)
                h = chebfun2(r.evalH(t));
            end
        end
        
        function f = get.p(r)
            f = @eval;
            function p = eval(t)
                p = chebfun2(r.evalP(t));
            end
        end
        
        %% Fluid volume computation
        function v = volume(r,t,H)
            [~,wx] = chebpts(r.n(1));  [~,wy] = chebpts(r.n(2));
            v = zeros(size(t));
            for k = 1:length(t)
                
                if nargin==2
                    H = r.evalH(t(k));
                end
                
                v(k) = volume_vals(r,t(k),H,wx,wy);
            end
        end
        
        %% Actual minimum value (from Cheb interpolation)
        function hm = minH(r,t)
            hm = zeros(length(t),1);
            for j = 1:length(t)
                hm(j) = min2( chebfun2(r.evalH(t(j))) );
            end
        end
        
        %% 2D plot functions
        function [XE,YE,H,P] = plotdata(r,t)
            F = chebfun2(r.evalH(t)');
            [X,Y] = ndgrid(linspace(-1,1,80),linspace(-1,1,60));
            H = F(X,Y);
            [XE,YE] = r.disc.map(t,X,Y);
            
            if nargout > 3
                F = chebfun2(r.evalP(t)');
                P = F(X,Y);
            end
        end
        
        function plot2d(r,t,var)
            if nargin < 3 || isequal(lower(var),'h')
                [XE,YE,U] = r.plotdata(t);
                cax = [0 r.Hmax];
            else
                [XE,YE,~,U] = r.plotdata(t);
                cax = [r.Pmin r.Pmax];
            end
            pcolor(XE,YE,U), shading interp
            axis equal, axis([-3 3 -0.8 0.8])
            %colormap redblue
            caxis(cax)
            xlabel('x'), ylabel('y')
            ti = sprintf('t = %.3f, S = %.2e, A = %.2e',t,r.pS,r.pA);
            title(ti)
            colorbar
        end
        
        function plot(r,t,var)
            if nargin < 3
                var = 'h';
            end
            [x,y] = r.shape(t);
            u = r.(lower(var));
            plot(x,y,u(t))
            %set(gca,'dataaspectratio',[1 1 6])
        end
        
        %% Save result
        % Saving the entire ODE solution structure is often prohibitively large.
        function savedata(r,fname,numt)
            if nargin < 3, numt = 250; end
            
            t = r.times(numt);
            
            dim = [ size(r.disc.square.grid.X) length(t) ];
            H = zeros(dim,'single');
            P = zeros(dim,'single');
            X = zeros(dim,'single');
            Y = zeros(dim,'single');
            for j = 1:length(t)
                H(:,:,j) = r.evalH(t(j));
                [X(:,:,j),Y(:,:,j)] = r.grid(t(j));
                P(:,:,j) = r.evalP(t(j));
            end
            
            volume = r.volume(t);
            r.solution = [];  result = struct(r);
            save(fname,'t','volume','X','Y','H','P','result');
        end
        
        %% Animation functions
        function animate(r,fname,t)
            figure
            r.plot2d(0)
            shg, hold on
            
            if nargin > 1
                vw = VideoWriter(fname,'MPEG-4');
                vw.FrameRate = 60;
                open(vw);
            end
            
            if nargin < 3
                t = r.times(250);
            end
            
            for j = 1:length(t)
                cla
                r.plot2d(t(j))
                if nargin > 1, writeVideo(vw,getframe(gcf));  end
                drawnow %pause(0.01)
            end
            
            if nargin > 1, close(vw); end
        end
        
        function slice(r,xx,fname)
            figure
            phimax = double(min(r.Hmax,16));
            if nargin < 2
                error('Provide a vector of x values in the square')
            end
            if nargin > 2
                vw = VideoWriter(fname,'MPEG-4');
                open(vw);
            end
            [~,yemax] = r.disc.eye.map.xy(0,1);  % for axis scaling
            t = r.times(250);
            hold on
            txt1 = text(.45,phimax*.92,'','fontsize',14,'hor','r');
            txt2 = text(.45,phimax*.82,'','fontsize',14,'hor','r');
            co = get(gca,'colororder');
            for j = 1:length(t)
                fun = chebfun(r.evalH(t(j)));  % chebfun in x at every y
                maxslope = -inf;  minval = inf;
                for kk = 1:length(xx)
                    % find where the x=const line ends up in the eye
                    [xl,yl] = r.disc.map(t(j),[xx(kk) xx(kk)],[-1 1]);
                    % make a chebfun out of values at the y nodes
                    u = chebfun( fun(xx(kk))', yl );
                    han(kk) = plot(u,'linewidth',2,'color',co(kk,:));
                    newmax = max( max(diff(u)),-min(diff(u)) );
                    maxslope = max(maxslope,newmax);
                    minval = min( minval, min(u) );
                end
                title(['t = ',num2str(t(j),'%.3f')])
                axis([-yemax yemax 0 phimax])
                xlabel('y'), ylabel('h'), shg
                txt1.String = sprintf('max slope: %.1f',maxslope);
                txt2.String = sprintf('min value: %.1e',minval);
                if nargin>2, writeVideo(vw,getframe(gcf));  end
                drawnow %pause(0.01)
                if j < length(t), delete(han), end
            end
            
            if nargin > 2, close(vw); end
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
            
            Motion = (r.disc.strip.dydt(t)'/r.disc.strip.map.deriv.yinv(t)) .* imag(gradHs);
            H_t = -DivQ - Motion;
            P_t = P + r.pS*DivH + r.pA*H.^-3;
            
            % Flux condition on left side (drainage)
            xleft = XS(1,:);  yleft = YS(1,:);
            absfp = r.disc.eye.map.absderiv(xleft,yleft);
            Qdotn = -Psi(1,:).*real(gradPs(1,:))./absfp;
            Qout = r.Qleft(t,yc);
            P_t(1,:) = Qdotn - Qout(:)';  % east
            
            % No flux on two sides
            xb = XS(end,:);  yb = YS(end,:);
            absfp = r.disc.eye.map.absderiv(xb,yb);
            Qdotn = Psi(end,:).*real(gradPs(end,:))./absfp;
            P_t(r.n(1),:) = Qdotn;   % west

            xb = XS(:,1);  yb = YS(:,1);
            absfp = r.disc.eye.map.absderiv(xb,yb);
            Qdotn = -Psi(:,1).*imag(gradPs(:,1))./absfp;
            P_t(:,1) = Qdotn;                % south
            
            % Flux on moving side
            xtop = XS(:,end);  ytop = YS(:,end);
            absfp = r.disc.eye.map.absderiv(xtop,ytop);
            Qdotn = Psi(:,r.n(2)).*imag(gradPs(:,r.n(2)))./absfp;
            Qin = r.Qtop(t,xc);
            P_t(:,r.n(2)) = Qdotn - Qin + r.drho_dt(t).*absfp.*(H(:,r.n(2))-r.h_e/2);
            
            u_t = r.pack(H_t,P_t); 
           % fprintf('\b\b\b\b\b\b%.4f',t)

        end
        
        function J = jac(r,t,u)
            [H,P] = r.unpack(u);
            [XS,YS,factor] = r.disc.strip.grid(t);
            Factor = r.disc.Diag(factor);
            
            [~,~,JH,Jgrads] = r.grad(t,H,factor);
            [gradP,~,JP] = r.grad(t,P,factor);
            
            Psi = H.^3/3;
            dPsi_dh = H.^2;
            
            %Q = -Psi.*gradP;
            DiagPsi = r.disc.Diag(Psi);
            J1 = [-r.disc.Diag(dPsi_dh.*gradP), -DiagPsi*JP ];
            
            %DivQ = div(t,Q,factor);
            J1 = real( Factor*conj(Jgrads)*J1 );
            %Motion = bsxfun(@times, dyc_dt(t)'./dyc_dys(t), imag(gradHs) );
            v = r.disc.strip.dydt(t)/r.disc.strip.map.deriv.yinv(t);
            JMo = kron( spdiags(v,0,r.n(2),r.n(2))*r.disc.strip.dm.y(t), r.disc.eye.x );
            %H_t = -DivQ - Motion;
            J1 = -J1 - [JMo,0*JMo];
            
            %DivH = div(t,gradH,factor);
            J2 = [real( Factor*conj(Jgrads)*JH ),0*JH];
            %P_t = P + pS*DivH + pA*H.^-3;
            J2 = [ -3*r.pA*r.disc.Diag(H.^(-4)), r.disc.eye.all ] + r.pS*J2;
            
            % Flux conditions
            %P_t([1 n(1)],:) = real(gradPs([1 n(1)],:));  % east and west
            %P_t(:,1) = imag(gradPs(:,1));                % south
            Jxs = real(Jgrads);  Jys = imag(Jgrads);
            bdy = r.disc.boundary;
            
            Jn = DiagPsi*Jxs;
            xleft = XS(1,:);  yleft = YS(1,:);
            absfp = r.disc.eye.map.absderiv(xleft,yleft);
            J2(bdy.E,:) = [0*Jxs(bdy.E,:) -diag(1./absfp)*Jn(bdy.E,:)];
            
            xb = XS(end,:);  yb = YS(end,:);
            absfp = r.disc.eye.map.absderiv(xb,yb);
            J2(bdy.W,:) = [0*Jxs(bdy.W,:) diag(1./absfp)*Jn(bdy.W,:)];
            
            Jn = DiagPsi*Jys;
            xb = XS(:,1);  yb = YS(:,1);
            absfp = r.disc.eye.map.absderiv(xb,yb);
            J2(bdy.S,:) = [0*Jys(bdy.S,:) -diag(1./absfp)*Jn(bdy.S,:)];
            
            % Moving side
            xtop = XS(:,end);  ytop = YS(:,end);
            absfp = r.disc.eye.map.absderiv(xtop,ytop);
            %Qdotn = Psi(:,n(2)).*imag(gradPs(:,n(2)))./absfp;
            %P_t(:,n(2)) = Qdotn + drho_dt(t).*absfp.*(H(:,n(2))-He);
            J2(bdy.N,:) = [diag(absfp)*r.drho_dt(t)*r.disc.eye.all(bdy.N,:),...
                diag(1./absfp)*Jn(bdy.N,:) ];
            %Qdotn = Pphi(:,n(2)).*Phi_ys(:,n(2))./amp;
            %Phi_t(:,n(2)) = Qdotn + drho_dt(t).*amp.*Phi(:,n(2));
            
            J=full([J1;J2]);
            
            % Boundary h values do not appear
            idx = find(r.disc.boundary.loc_outer);
            J(idx,:) = [ ];    J(:,idx) = [ ];
            
        end
        
        function [H,P,du] = initial(r)
            g = r.disc.square.grid;
            IC = r.initcond;
            du = [];
            if isa(IC,'struct')
                H = pack(r,IC.H(g.X,g.Y),IC.P(g.X,g.Y));
                P = pack(r,IC.dH(g.X,g.Y),IC.dP(g.X,g.Y));
                du = 1;
            elseif isnumeric(IC)
                H = IC(:,1);
                P = IC(:,2);
                du = 1;
            else
                switch r.initcond
                    case 'flat'
                        H = r.boundaryH*ones(r.n);
                    case 'laplace'
                        bidx = r.disc.boundary.all;
                        D = r.disc.square.dm;   I = r.disc.eye;
                        L = kron(D.y^2,I.x) + kron(I.y,D.x^2);
                        L(bidx,:) = 0;
                        nb = sum(bidx(:));
                        L(bidx,bidx) = eye(nb);
                        b = 10*(g.X.^4+g.Y.^4);  b = b(:);
                        b(bidx) = 1;
                        H = reshape(L\b,r.n);
                        
                        % affine xform to get the right initial volume and BC
                        A = [1 H(1); r.volume_vals(0,H.^0) r.volume_vals(0,H)];
                        c = A \ [r.boundaryH;r.initvolume];
                        H = c(1) + c(2)*H;
                        
                        %hmin = min(H(:));
                        %H = r.boundaryH*((H-hmin)/(1-hmin)*.95 + .05);
                end
                
                [~,~,factor] = r.disc.strip.grid(0);
                gradH = grad(r,0,H,factor);
                LapH = div(r,0,gradH,factor);
                P = -r.pA./H.^3 - r.pS*LapH;
                
            end
            
        end
        
        function v = volume_fun_subgrid(r,t,U,xc,yc,wx,wy)
                        
            
            
            A = (r.disc.strip.map.jaccs_xx(xc).*(r.rho(t)+1)/2) .* U;
            [XS,YS] = r.disc.strip.grid(t,xc,yc);
            % Jacobian for strip to eye map
            A = A.*r.disc.eye.map.absderiv(XS,YS).^2;
            v = (wx*A*wy');
            
            
            v = (wx*A*wy');
            
        end
        
    end
    
    methods (Access=private)
        function r = setup_discretization(r)
            
            % Discretization in [-1,1]^2
            xc = chebpts(r.n(1),r.domain(1,:));  Dxc = diffmat(r.n(1),1,r.domain(1,:));
            yc = chebpts(r.n(2),r.domain(2,:));  Dyc = diffmat(r.n(2),1,r.domain(2,:));
            [XC,YC] = ndgrid(xc,yc);
            Ix = speye(r.n(1));  Iy = speye(r.n(2));
            I = speye(prod(r.n));
            
            % Transformation from square (xc,yc) to strip (xs,ys)
            alpha = 1.4;  gamma = 7;
            xc2xs = @(xc) gamma*xc./(alpha^2-xc.^2);
            yc2ys = @(t,yc) (yc+1).*(1+r.rho(t))/2-1;  % map from [-1,1] to [-1,ymax(t)]
            
            jaccs_xx = @(x)  gamma*(alpha^2+x.^2)./(alpha.^2-x.^2).^2;
            discr.jaccs_xx = jaccs_xx;
            dxs_dxc = jaccs_xx(xc);
            Dxs = diag(1./dxs_dxc)*Dxc;
            jacsc_yy = @(t) 2./(1+r.rho(t));
            dyc_dys = @(t) 2./(1+r.rho(t));
            Dys = @(t) Dyc*dyc_dys(t);  % leaves out the d/dt term for below
            dyc_dt = @(t) -r.drho_dt(t)*(yc+1)/(r.rho(t)+1);
            
            function [XS,YS,factor] = stripgrid(t,xc_o,yc_o)
                
                if nargin==1
                    [XS,YS] = ndgrid(xc2xs(xc),yc2ys(t,yc));
                else
                    [XS,YS] = ndgrid(xc2xs(xc_o),yc2ys(t,yc_o));
                end
                
                if nargout > 2
                    factor = 1./dze_dzs(XS,YS);
                end
                
            end
            
            
            
            % Transformation from strip (xs,ys) to eye (xe,ye)
            mapc = 3.84;
            s2e = @(xs,ys) deal(3.0*real(tanh((xs+1i*ys)/mapc)),3.0*imag(tanh((xs+1i*ys)/mapc)));
            dze_dzs = @(xs,ys) 3.0*sech((xs+1i*ys)/mapc).^2 / mapc;
            abs_dze_dzs = @(xs,ys) 3.0*(2/mapc)*(cosh(2*xs/mapc)+cos(2*ys/mapc)).^(-1);
            
            c2e = @(t,xc,yc) s2e( xc2xs(xc), yc2ys(t,yc) );  % composite c->e
            
            % Identify the boundaries in the grid
            east = false(r.n);   east(1,:) = true;
            west = false(r.n);   west(r.n(1),:) = true;
            south = false(r.n);  south(:,1) = true;
            north = false(r.n);  north(:,r.n(2)) = true;
            
            discr.eye = struct('x',Ix,'y',Iy,'all',I);
            discr.vec = @(U) U(:);  discr.unvec = @(u) reshape(u,r.n);
            discr.Diag = @(U) spdiags(discr.vec(U),0,prod(r.n),prod(r.n));
            discr.num = struct('h',prod(r.n-2),'p',prod(r.n));
            discr.square.points = struct('x',xc,'y',yc);
            discr.square.grid = struct('X',XC,'Y',YC);
            discr.square.dm = struct('x',Dxc,'y',Dyc);
            discr.strip.map = struct('x',xc2xs,'y',yc2ys);
            discr.strip.map.deriv = struct('x',dxs_dxc,'yinv',dyc_dys);
            discr.strip.map.jaccs_xx = jaccs_xx;
            discr.strip.map.jacsc_yy = jacsc_yy;
            discr.strip.points = struct('x',xc2xs(xc),'y',@(t) yc2ys(t,yc));
            discr.strip.grid = @stripgrid;
            discr.strip.dydt = dyc_dt;
            discr.strip.dm = struct('x',Dxs,'y',Dys);
            discr.eye.map = struct('xy',s2e,'deriv',dze_dzs,'absderiv',abs_dze_dzs);
            discr.map = c2e;
            discr.boundary = struct('N',north,'S',south,'E',east,'W',west);
            discr.boundary.all = east | west | north | south;
            discr.boundary.loc_outer = discr.boundary.all;
            
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
        

        
        function st = ode_status(r,t,y,flag)
            persistent wb
            if isequal(flag,'init')
                wb = waitbar(0,'Opening phase','name','Progress');
            elseif isequal(flag,'done')
                delete(wb)
            else
                tt = rem(t,r.period);
                if tt < 0.19
                    waitbar(tt/.19,wb,'Opening phase');
                elseif tt < 5.176
                    waitbar((tt-0.19)/(5.176-0.19),wb,'Open phase');
                else
                    waitbar((tt-5.176)/(r.period-5.176),wb,'Closing phase');
                end
            end
            st = 0;
        end
        
        function [Qtop,Qleft] = fluxfuns(r)
            
            period = r.period;
            
            %% domain mapping functions
            dxs_dxc = r.disc.strip.map.jaccs_xx;
            xc2xs = r.disc.strip.map.x;
            yc2ys = r.disc.strip.map.y;
            dyc_dys = r.disc.strip.map.deriv.yinv;
            abs_dze_dzs = r.disc.eye.map.absderiv;
                  
            %% shape/windowing functions
            humpfun = @(x,center,hw) exp( -log(2)*((x-center)/hw).^2 );
            function v = stepfun(x,start,stop)
                c = (stop+start)/2;
                d = (stop-start)/2;
                v = (1 + tanh(4*(x-c)/d))/2;
            end
            function w = windowfun(x,start,thick,stop)
                w = stepfun(x,start,start+thick) - stepfun(x,stop-thick,stop);
            end
            
            %% Pointwise flux evaluations
             function Q = qtop(t,xc)
                ts = 0.02;   % start time for the window function
                tth = 0.1;   % tanh thickness
                tfrac = 0.3;  % fraction of blink cycle time
                xcen = 0.5;  % center of hump along the top
                xwid = 0.25; % width of the hump

                Q = windowfun(t,ts,tth,tfrac*period).*humpfun(xc,xcen,xwid);
                
%                Jcs = dxs_dxc(xc);
%                Jse = abs_dze_dzs(xc2xs(xc),yc2ys(t,1));
%                Q = Q.*Jcs.*Jse;
             end
            
            function Q = qleft(t,yc)
                ts = 0.05;   % starting time
                tth = 0.1;   % thickness
                tfrac = 0.5;  % fraction of blink cycle time
                % space values
                % need decay by ends of domain for Gaussian (hump)
                ys = -0.9;   % starting pos
                yth = 0.25;  % thickness
                ye = 0.9;    % end position

                Q = windowfun(t,ts,tth,tfrac*period).*windowfun(yc,ys,yth,ye);
                
%                Jcs = 1./dyc_dys(t);
%                Jse = abs_dze_dzs(xc2xs(-1),yc2ys(t,yc));
%                Q = Q.*Jcs.*Jse;
            end
                                   
            %% Integrate flux along one edge, at one time
            function S = flux_top(t,qc)
                [xc,wght] = chebpts(80);
                Jcs = dxs_dxc(xc);
                xs = xc2xs(xc);
                Jse = abs_dze_dzs(xs,yc2ys(t,1));
                S = sum( qc(t,xc).*Jcs.*Jse.*wght' );
%                S = sum( qc(t,xc).*wght' );
            end
            
            function S = flux_bot(t,qc)
                [xc,wght] = chebpts(80);
                Jcs = dxs_dxc(xc);
                xs = xc2xs(xc);
                Jse = abs_dze_dzs(xs,yc2ys(t,-1));
                S = sum( qc(t,xc).*Jcs.*Jse.*wght' );
            end
            
            function S = flux_left(t,qc)
                [yc,wght] = chebpts(80);
                 Jcs = 1/dyc_dys(t);
                 ys = yc2ys(t,yc);
                 Jse = abs_dze_dzs(xc2xs(-1),ys);
                 S = sum( qc(t,yc).*Jcs.*Jse.*wght' );
%                S = sum( qc(t,yc).*wght' );
            end
            
            %% Integrate flux along an edge, over one blink cycle. 
            function Q = totalflux(flux,q)
                [node,wght] = chebpts(80);
                t = (node+1)/2*period;
                S = zeros(size(t));
                for ii = 1:length(t)
                    S(ii) = flux(t(ii),q);
                end
                Q = period/2*sum( wght'.*S );
            end
            
            %% Find amplitudes to meet supply/drainage target
            Ain = r.supplyvolume / totalflux(@flux_top,@qtop);
            Qtop = @(t,xc) Ain*qtop(t,xc);
            
            Aout = r.drainvolume / totalflux(@flux_left,@qleft);
            Qleft = @(t,xc) -Aout*qleft(t,xc);

        end
        
    end
    
    
    
    
end