classdef blinkmulti
	
	properties
		region  % array of blinkone objects
		H       % PUC tree
		P       % PUC tree
		
		model   % parameters of the model
		map     % maps between the domains
		
		tspan = [0 5.258]
		
		method   % "NKS" or "SNK"
		odetol = 1e-4
		solvetime
		solution
		finalstate
		initstate
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
		
		function b = blinkmulti(model,space,time)
			% Require Chebfun
			assert(exist('chebpts','file')>0,'Add Chebfun to path');
			
			%*** model parameters
			pp = inputParser;
			pp.KeepUnmatched = true;
			pp.addParameter('A',NaN);
 			pp.addParameter('S',NaN);
 			pp.addParameter('h_slideover',0);
 			pp.addParameter('h_boundary',13);
 			pp.addParameter('percentclosed',0);
 			pp.addParameter('drainvolume',0);
			pp.addParameter('supplyvolume',0);
			parse(pp,model);
			b.model = pp.Results;
			
			%*** spatial parameters and setup
			ps = inputParser;
			ps.KeepUnmatched = true;
			ps.addParameter('degree',[10 10]);
			ps.addParameter('coarsedegree',[5 5]);
			ps.addParameter('tol',1e-5);
			ps.addParameter('splitdim',2);
			parse(ps,space);
			ps = ps.Results;
			
			b.percentClosed = b.model.percentclosed;  % trigger setting rho
			b = coord_maps(b);
					
			% Create tree by repeated splits
			T = make_tree(b,ps.degree,ps.coarsedegree,ps.tol,ps.splitdim);

			% solution data structures
			b.H = PUchebfun(T);
			setInterpMatrices(b.H,false);
			b.P = copy(b.H);
			
			% set up subdomain objects
			Hleaf = b.H.leafArray;
			for i = 1:length(Hleaf)
				b.region{i} = blinkone(b.model,Hleaf{i},b.map);
			end			

			%*** temporal parameters
			pt = inputParser;
			pt.KeepUnmatched = true;
			pt.addParameter('tol',1e-4);
			pt.addParameter('tspan',[0 5.258]);  % hard-coded period
			pt.addParameter('method',"NKS");
			pt.addParameter('initstate',[]);
			parse(pt,time);
			pt = pt.Results;
			
			b.method = pt.method;
			b.odetol = pt.tol;
			b.initstate = pt.initstate;
			b.tspan = pt.tspan;
				
			% use constant initial state if none given
			if isempty(b.initstate)
				b.initstate.H = chebfun2(b.model.h_boundary);
				b.initstate.P = chebfun2(0);
				b.initstate.dH = chebfun2(0);
				b.initstate.dP = chebfun2(0);
			end
			
			b.the_period = b.region{1}.period;
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
		
		function l = length(b)
			l = length(b.region);
		end
		
		function [b,M,u0,du0] = initialize(b)
			% mass matrices
			m = length(b);
			M = cell(m,1);
			for i = 1:m
				r = b.region{i};
				MH = ones(r.disc.num.h,1);
				MP = zeros(r.disc.num.p,1);
				N = r.disc.num.h + r.disc.num.p;
				M{i} = spdiags([MH;MP],0,N,N);
			end
			
			state = b.initstate;
			
			% initial slope
			dH = copy(b.H);
			dP = copy(b.P);
			dH.sample(state.dH);
			dP.sample(state.dP);
			dH.pack();  % get the right unknowns this way
			du0 = [ dH.Getvalues(); dP.Getvalues() ];
			
			% initial value
			b.H.sample(state.H);
			b.P.sample(state.P);
			b.H.pack();
			u0 = [ b.H.Getvalues(); b.P.Getvalues() ];
		end
		
		function b = solve(b,tspan)
			[b,M,u0,du0] = initialize(b);
			
			%% options
			opt = odeset('mass',M,...
				'reltol',b.odetol,'abstol',b.odetol,...
				'initialslope',du0,...
				'outputfcn',@(t,y,f) ode_status(b,t,y,f)   );
			if nargin > 1
				b.tspan = tspan;
			end

			%% solve			
			fprintf('\n\nt = 0.0000')
			tic
			sol = ASode15s(b.method=="SNK",...
				b.region,...
				b.tspan,...
				u0,...
				{b.H,b.P},...
				1,opt);
			b.solvetime = toc;
			fprintf('\n\n')
			
			%% extract final state (for solution continuation)
			b.solution = sol;
			b.finalstate.H = evalH(b,sol.x(end));
			b.finalstate.P = evalP(b,sol.x(end));
			
			% slope is the tricky bit, due to boundary
			[~,dufinal] = deval(sol,sol.x(end));
			dH = copy(b.H);
			for i = 1:length(b)
				dH.leafArray{i}.SetBoundaryValues(0);
			end
			dP = copy(b.P);
			b.finalstate.dH = evalsol(b,dufinal(1:length(b.H)),dH);
			b.finalstate.dP = evalsol(b,dufinal(length(b.H)+1:end),dP);
			
			%% save memory
			if b.odetol > 1e-7
				sol.y = single(sol.y);
				sol.idata.dif3d = single(sol.idata.dif3d);
			end
			b.solution = sol;
			
			% may not have finished, and rounding to single can "change"
			% the final time
			b.tspan = [sol.x(1) sol.x(end)];
			
		end
		
		function v = get.Hmin(b)
			% Minimum value over the collocation values and all time.
			r = b.region{1};
			v = min(r.boundaryH,min(min( b.solution.y(1:r.disc.num.h,:) )));
		end
		
		function v = get.Hmax(b)
			% Maximum value over the collocation values and all time.
			r = b.region{1};
			v = max(r.boundaryH,max(max( b.solution.y(1:r.disc.num.h,:) )));
		end
		
		function v = get.Pmin(b)
			% Minimum value over the collocation values and all time.
			r = b.region{1};
			v = min(min( b.solution.y(1+r.disc.num.h:end,:) ));
		end
		
		function v = get.Pmax(b)
			% Maximum value over the collocation values and all time.
			r = b.region{1};
			v = max(max( b.solution.y(1+r.disc.num.h:end,:) ));
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
		
		function [X,Y] = grid(b,t)
			% grid of points in the eye domain
			[X,Y] = meshgrid(chebpts(80));
			[X,Y] = b.region{1}.disc.map(t,X,Y);
		end
		
		function [x,y] = shape(b,t)
			% chebfuns for the x and y coordinates in the eye
			[X,Y] = b.grid(t);
			x = chebfun2(X);
			y = chebfun2(Y);
		end
				
		function H = evalH(b,t)
			u = deval(b.solution,t);
			H = evalsol(b,u(1:length(b.H)),b.H);
		end
		
		function P = evalP(b,t)
			u = deval(b.solution,t);
			P = evalsol(b,u(length(b.H)+1:end),b.P);
		end
		
		function f = get.h(b)
			f = @(t) evalH(b,t);
		end
		
		function f = get.p(b)
			f = @(t) evalP(b,t);
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
			F = r.evalH(t);
			[X,Y] = ndgrid(linspace(-1,1,90),linspace(-1,1,80));
			H = F(X,Y);
			[XE,YE] = r.region{1}.disc.map(t,X,Y);
			
			if nargout > 3
				F = r.evalP(t);
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
			ti = sprintf('t = %.3f, S = %.2e, A = %.2e',t,r.region{1}.pS,r.region{1}.pA);
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
		
		
		%% Helper for the H and P evaluation functions
		function v = evalsol(b,u,V)
			sample(V,u);			
			n = 40; 
			x = chebpts(n);
			v = chebfun2( evalfGrid(V,{x,x})' );
			done = false;
			while ~done
				n = ceil(1.4*n);
				vold = v;
				x = chebpts(n);
				v = chebfun2( evalfGrid(V,{x,x})' );
				if n > 300
					warning("can't resolve solution with a chebfun")
					break
				end
				done = rank(v)==rank(vold) || norm(v-vold) > b.odetol/4*norm(v);
			end
		end
		
		%% Functions for obtaining the solution.
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
		
		function T = make_tree(b,degree,cdegree,tol,splitdim)
			param = struct('degs',degree,...
				'cdegs',cdegree,...
				'split_flag',[true, true],...
				'tol',tol,...
				'domain',[-1 1;-1 1]  );
			T = ChebPatch(param);
			expr = "T";
			for k = 1:length(splitdim)
				s = splitdim(k);
				for j = 1:length(expr)
					eval(expr(j) + " = split(" + expr(j) + "," + s + ");")
				end
				expr = [ expr+".children{1}", expr+".children{2}" ];
				expr = expr(:);
			end
			clean(T);
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
		
		function b = coord_maps(b)
			% From square "c" to strip "s"
			alpha = 1.4;  gamma = 7;
            b.map.strip.x = @(xc) gamma*xc./(alpha^2-xc.^2);
            b.map.strip.y = @(t,yc) (yc+1).*(1+b.rho(t))/2-1;  % map from [-1,1] to [-1,ymax(t)]
			b.map.strip.dxdx = @(x)  gamma*(alpha^2+x.^2)./(alpha.^2-x.^2).^2;  % dx_s/dx_c
            b.map.strip.dydyinv = @(t) 2./(1+b.rho(t));       % dy_c/dy_s
            b.map.strip.dydt = @(t,yc) -b.drho_dt(t)*(yc+1)/(b.rho(t)+1);   % dy_s/dt
			
			% From strip "s" to eye "e"
            mapc = 3.84;
            b.map.eye.xy = @(xs,ys) deal(3.0*real(tanh((xs+1i*ys)/mapc)),3.0*imag(tanh((xs+1i*ys)/mapc)));
            b.map.eye.dzdz = @(xs,ys) 3.0*sech((xs+1i*ys)/mapc).^2 / mapc;
            b.map.eye.absdzdz = @(xs,ys) 3.0*(2/mapc)*(cosh(2*xs/mapc)+cos(2*ys/mapc)).^(-1);
            
			% Composite from "c" to "e"
            b.map.both = @(t,xc,yc) b.map.eye.xy( b.map.strip.x(xc), b.map.strip.y(t,yc) );  % composite c->e
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