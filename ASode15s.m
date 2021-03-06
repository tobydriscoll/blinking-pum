function varargout = ASode15s(solveoptions,ode,tspan,y0,PUApprox,interface_scale,options,varargin)

DEBUG = 3;
MAXITER = 6;

solver_name = 'ode15s';

if nargin < 7
    options = [];
end

use_SNK = solveoptions.method == "SNK";
use_parallel = solveoptions.use_parallel;

% Stats
nsteps   = 0;
nfailed  = 0;
nfevals  = 0;
npds     = 0;
ndecomps = 0;
nsolves  = 0;

% Output
FcnHandlesUsed  = true;   % isa(ode,'function_handle');
output_sol = (FcnHandlesUsed && (nargout==1));      % sol = odeXX(...)
output_ty  = (~output_sol && (nargout > 0));  % [t,y,...] = odeXX(...)
% There might be no output requested...

sol = []; kvec = []; dif3d = [];
if output_sol
    sol.solver = solver_name;
    sol.extdata.odefun = ode;
    sol.extdata.options = options;
    sol.extdata.varargin = varargin;
end

% Handle solver arguments
[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, odeFcn, ...
    options, threshold, rtol, normcontrol, normy, hmax, htry, htspan] = ...
    odearguments(FcnHandlesUsed, solver_name, ode, tspan, y0, options, varargin);
nfevals = nfevals + 1;
normcontrol = true;
normy = norm(y0)/length(y0);
%This won't work with row scale! needs to be local to each patch
one2neq = (1:neq);

% Handle the output
if nargout > 0
    outputFcn = odeget(options,'OutputFcn',[],'fast');
else
    outputFcn = odeget(options,'OutputFcn',@odeplot,'fast');
end
outputArgs = {};
if isempty(outputFcn)
    haveOutputFcn = false;
else
    haveOutputFcn = true;
    outputs = odeget(options,'OutputSel',1:neq,'fast');
    if isa(outputFcn,'function_handle')
        % With MATLAB 6 syntax pass additional input arguments to outputFcn.
        outputArgs = varargin;
    end
end
refine = max(1,odeget(options,'Refine',1,'fast'));
if ntspan > 2
    outputAt = 'RequestedPoints';         % output only at tspan points
elseif refine <= 1
    outputAt = 'SolverSteps';             % computed points, no refinement
else
    outputAt = 'RefinedSteps';            % computed points, with refinement
    S = (1:refine-1) / refine;
end
printstats = strcmp(odeget(options,'Stats','off','fast'),'on');

% Handle the event function NO!
[haveEventFcn,eventFcn,eventArgs,valt,teout,yeout,ieout] = ...
    odeevents(FcnHandlesUsed,odeFcn,t0,y0,options,varargin);
%haveEventFcn = false;
%teout = [];
%yeout = [];
%ieout = [];

% Handle the mass matrix NO! Mtype = 1 (constant). Assume cell array of mass matrices
% is passed through (one for each patch).
%[Mtype, Mt, Mfun, Margs, dMoptions] = odemass(FcnHandlesUsed,odeFcn,t0,y0,...
%                                             options,varargin);

Mt = odeget(options,'Mass',[],'fast');
Mtype = 1;

% Non-negative solution components NO! We don't worry about this
%
idxNonNegative = odeget(options,'NonNegative',[],'fast');

% HEY! I dont think this is needed
% Handle the Jacobian
[Jconstant,Jac,Jargs,Joptions] = ...
    odejacobian(FcnHandlesUsed,odeFcn,t0,y0,options,varargin);
Janalytic = isempty(Joptions);

t = t0;
y = y0;

yp0_OK = true;
DAE = true;
RowScale = [];
if Mtype > 0
    %   nz = nnz(Mt);
    %   if nz == 0
    %     error(message('MATLAB:ode15s:MassMatrixAllZero'))
    %   end
    
    Msingular = 'yes'; DAE = true;
    %   Msingular = odeget(options,'MassSingular','maybe','fast');
    %   switch Msingular
    %     case 'no',     DAE = false;
    %     case 'yes',    DAE = true;
    %     case 'maybe',  DAE = (eps*nz*condest(Mt) > 1);
    %   end
    
    if DAE
        yp0 = odeget(options,'InitialSlope',[],'fast');
        %if isempty(yp0)
        % yp0_OK = false;
        % yp0 = zeros(neq,1);
        %else
        yp0 = yp0(:);
        % if length(yp0) ~= neq
        %   error(message('MATLAB:ode15s:YoYPoLengthMismatch'));
        % end
        % NO! we assume things are fine
        % Test if (y0,yp0) are consistent enough to accept.
        %yp0_OK = (norm(Mt*yp0 - f0) <= 1e-3*rtol*max(norm(Mt*yp0),norm(f0)));
        %end
        % NO! assume everything is fine
        %     if ~yp0_OK           % Must compute ICs, so classify them.
        %       if Mtype >= 3  % state dependent
        %         ICtype = 3;
        %       else  % M, M(t)
        %         % Test for a diagonal mass matrix.
        %         [r,c] = find(Mt);
        %         if isequal(r,c)   % diagonal
        %           ICtype = 1;
        %         elseif ~issparse(Mt) % not diagonal but full
        %           ICtype = 2;
        %         else  % sparse, not diagonal
        %           ICtype = 3;
        %         end
        %       end
        %     end
    end
end
Mcurrent = true;
Mtnew = Mt;

% if not set via 'options', initialize constant Jacobian here
% Maybe provide function for jacobian? IDK. needs to be patch by patch
% if Jconstant
%   if isempty(Jac) % use odenumjac
%     [Jac,Joptions.fac,nF] = odenumjac(odeFcn, {t0,y0,odeArgs{:}}, f0, Joptions);
%     nfevals = nfevals + nF;
%     npds = npds + 1;
%   elseif ~isa(Jac,'numeric')  % not been set via 'options'
%     Jac = feval(Jac,t0,y0,Jargs{:}); % replace by its value
%     npds = npds + 1;
%   end
% end

maxk = odeget(options,'MaxOrder',5,'fast');
bdf = strcmp(odeget(options,'BDF','off','fast'),'on');

% Initialize method parameters.
G = [1; 3/2; 11/6; 25/12; 137/60];
if bdf
    alpha = [0; 0; 0; 0; 0];
else
    alpha = [-37/200; -1/9; -0.0823; -0.0415; 0];
end
invGa = 1 ./ (G .* (1 - alpha));
erconst = alpha .* G + (1 ./ (2:6)');
difU = [ -1, -2, -3, -4,  -5;           % difU is its own inverse!
    0,  1,  3,  6,  10;
    0,  0, -1, -4, -10;
    0,  0,  0,  1,   5;
    0,  0,  0,  0,  -1 ];
maxK = 1:maxk;
[kJ,kI] = meshgrid(maxK,maxK);
difU = difU(maxK,maxK);

% Adjust the warnings.
warnoffId = { 'MATLAB:singularMatrix', 'MATLAB:nearlySingularMatrix'};
for i = 1:length(warnoffId)
    warnstat(i) = warning('query',warnoffId{i});
    warnoff(i) = warnstat(i);
    warnoff(i).state = 'off';
end


%yp0 = Masstimes(PUApprox,Mt,NKSresidual(t0,y0,1,PUApprox,ode,interface_scale));

if isempty(yp0)
	% should never be used
	[y,yp0] = GetInitialSlope(Mt,y0,zeros(size(y0)),t0,PUApprox,ode,rtol,interface_scale);
end

%yp0 = y0;

% NO! we assum yp0_ok = true
%
% Get the initial slope yp. For DAEs the default is to compute
% consistent initial conditions.
% if DAE && ~yp0_OK
%   if ICtype < 3
%     [y,yp,f0,dfdy,nFE,nPD,Jfac] = daeic12(odeFcn,odeArgs,t,ICtype,Mt,y,yp0,f0,...
%                                           rtol,Jconstant,Jac,Jargs,Joptions);
%   else
%     [y,yp,f0,dfdy,nFE,nPD,Jfac,dMfac] = daeic3(odeFcn,odeArgs,tspan,htry,Mtype,Mt,Mfun,...
%                                                Margs,dMoptions,y,yp0,f0,rtol,Jconstant,...
%                                                Jac,Jargs,Joptions);
%     if ~isempty(dMoptions)
%       dMoptions.fac = dMfac;
%     end
%   end
%   if ~isempty(Joptions)
%     Joptions.fac = Jfac;
%   end
%   nfevals = nfevals + nFE;
%   npds = npds + nPD;
%   if Mtype >= 3
%     Mt = feval(Mfun,t,y,Margs{:});
%     Mtnew = Mt;
%     Mcurrent = true;
%   end
%else

%   if Mtype == 0
%     yp = f0;
%   elseif DAE && yp0_OK
%     yp = yp0;
%
% % NO! DAE is true
% %
% %   else
% %     if issparse(Mt)
% %       [L,U,P,Q,R] = lu(Mt);
% %       yp = Q * (U \ (L \ (P * (R \ f0))));
% %     else
% %       [L,U,p] = lu(Mt,'vector');
% %       yp = U \ (L \ f0(p));
% %     end
% %     ndecomps = ndecomps + 1;
% %     nsolves = nsolves + 1;
%
%   end

% NO! change to match cell array of jacs.
%
%   if Jconstant
%     dfdy = Jac;
%   elseif Janalytic
%     dfdy = feval(Jac,t,y,Jargs{:});
%     npds = npds + 1;
%   else   % Joptions not empty
%     [dfdy,Joptions.fac,nF] = odenumjac(odeFcn, {t,y,odeArgs{:}}, f0, Joptions);
%     nfevals = nfevals + nF;
%     npds = npds + 1;
%   end
%end

yp = yp0;

%cell array of Jacobians
%dfdy = ComputeJac(Jac,num_sols,PUApprox,t,y);

Jcurrent = true;

% hmin is a small number such that t + hmin is clearly different from t in
% the working precision, but with this definition, it is 0 if t = 0.
hmin = 16*eps*abs(t);

if isempty(htry)
    % Compute an initial step size h using yp = y'(t).
    if normcontrol
        wt = max(normy,threshold);
        % yp_n = norm(Masstimes(PUApprox,Mt,yp));
        rh = 1.25 * (norm(yp) / wt) / sqrt(rtol);  % 1.25 = 1 / 0.8
    else
        wt = max(abs(y),threshold);
        % yp_n = norm(Masstimes(PUApprox,Mt,yp));
        rh = 1.25 * norm(yp./ wt,inf) / sqrt(rtol);
    end
    absh = min(hmax, htspan);
    if absh * rh > 1
        absh = 1 / rh;
    end
    absh = max(absh, hmin);
    
    %NO! we don't have to worry about this.
    %   if ~DAE
    %     % The error of BDF1 is 0.5*h^2*y''(t), so we can determine the optimal h.
    %     h = tdir * absh;
    %     tdel = (t + tdir*min(sqrt(eps)*max(abs(t),abs(t+h)),absh)) - t;
    %     f1 = feval(odeFcn,t+tdel,y,odeArgs{:});
    %     nfevals = nfevals + 1;
    %     dfdt = (f1 - f0) ./ tdel;
    %     DfDt = dfdt + dfdy*yp;
    %     if normcontrol
    %       if Mtype > 0
    %           if issparse(Mt)
    %               rh = 1.25 * sqrt(0.5 * (norm(U \ (L \ (P * (R \ DfDt)))) / wt) / rtol);
    %           else
    %               rh = 1.25 * sqrt(0.5 * (norm(U \ (L \ DfDt(p))) / wt) / rtol);
    %           end
    %       else
    %         rh = 1.25 * sqrt(0.5 * (norm(DfDt) / wt) / rtol);
    %       end
    %     else
    %       if Mtype > 0
    %         if issparse(Mt)
    %           rh = 1.25*sqrt(0.5*norm((Q * (U \ (L \ (P * (R \ DfDt))))) ./ wt,inf) / rtol);
    %         else
    %           rh = 1.25*sqrt(0.5*norm((U \ (L \ DfDt(p))) ./ wt,inf) / rtol);
    %         end
    %       else
    %         rh = 1.25 * sqrt(0.5 * norm( DfDt ./ wt,inf) / rtol);
    %       end
    %     end
    %     absh = min(hmax, htspan);
    %     if absh * rh > 1
    %       absh = 1 / rh;
    %     end
    %     absh = max(absh, hmin);
    %   end
else
    absh = min(hmax, max(hmin, htry));
end
h = tdir * absh;

% Initialize.
k = 1;                                  % start at order 1 with BDF1
K = 1;                                  % K = 1:k
klast = k;
abshlast = absh;

dif = zeros(neq,maxk+2);
dif(:,1) = h * yp;

hinvGak = h * invGa(k);
nconhk = 0;                             % steps taken with current h and k


%Miter = timeDiff(PUApprox,Mt,dfdy,hinvGak);
%Miter = Mt - hinvGak * dfdy;

% NO! we assume mass matrix is constant
%
% Account for strongly state-dependent mass matrix.
% if Mtype == 4
%   psi = dif(:,K) * (G(K) * invGa(k));
%   [dMpsidy,dMoptions.fac] = odenumjac(@odemxv, {Mfun,t,y,psi,Margs{:}}, Mt*psi, ...
%                                       dMoptions);
%   Miter = Miter + dMpsidy;
% end

% HEY! this needs to change.
% Use explicit scaling of the equations when solving DAEs.

%[RowScale,Miter] = RowScales(PUApprox,Miter,num_sols);
%[L,U,p] = LUarray(PUApprox,Miter);

% if DAE
%   RowScale = 1 ./ max(abs(Miter),[],2);
%   Miter = sparse(one2neq,one2neq,RowScale) * Miter;
% end
% if issparse(Miter)
%   [L,U,P,Q,R] = lu(Miter);
% else
%   [L,U,p] = lu(Miter,'vector');
% end


ndecomps = ndecomps + 1;
havrate = false;

% Allocate memory if we're generating output.
nout = 0;
tout = []; yout = [];
if nargout > 0
    if output_sol
        chunk = min(max(100,50*refine), refine+floor((2^11)/neq));
        tout = zeros(1,chunk);
        yout = zeros(neq,chunk);
        kvec = zeros(1,chunk);
        dif3d = zeros(neq,maxk+2,chunk);
    else
        if ntspan > 2                         % output only at tspan points
            tout = zeros(1,ntspan);
            yout = zeros(neq,ntspan);
        else                                  % alloc in chunks
            chunk = min(max(100,50*refine), refine+floor((2^13)/neq));
            tout = zeros(1,chunk);
            yout = zeros(neq,chunk);
        end
    end
    nout = 1;
    tout(nout) = t;
    yout(:,nout) = y;
end

% Initialize the output function.
if haveOutputFcn
    feval(outputFcn,[t tfinal],y(outputs),'init',outputArgs{:});
end

% THE MAIN LOOP
H_sol = y(1:length(PUApprox{1}));
PUApprox{1}.sample(H_sol);
%int_vol = BlinkVolume(ode,PUApprox{1},t);


done = false;
at_hmin = false;
while ~done
    
    hmin = 16*eps(t);
    absh = min(hmax, max(hmin, absh));
    if absh == hmin
        if at_hmin
            absh = abshlast;  % required by stepsize recovery
        end
        at_hmin = true;
    else
        at_hmin = false;
    end
    h = tdir * absh;
    
    % Stretch the step if within 10% of tfinal-t.
    if 1.1*absh >= abs(tfinal - t)
        h = tfinal - t;
        absh = abs(h);
        done = true;
    end
    
    if (absh ~= abshlast) || (k ~= klast)
        difRU = cumprod((kI - 1 - kJ*(absh/abshlast)) ./ kI) * difU;
        dif(:,K) = dif(:,K) * difRU(K,K);
        
        hinvGak = h * invGa(k);
        nconhk = 0;
        
        %cell array with local linear operators for time step integration
        %Miter = timeDiff(PUApprox,Mt,dfdy,hinvGak);
        
        ndecomps = ndecomps + 1;
        havrate = false;
    end
    
    min_iter = 1;
    inter_tol = 1e-10;
    res_tol = 1e-4;
    
    % LOOP FOR ADVANCING ONE STEP.
    nofailed = true;                      % no failed attempts
    while true                            % Evaluate the formula.
        
        gotynew = false;                    % is ynew evaluated yet?
        while ~gotynew
            
            % Compute the constant terms in the equation for ynew.
            psi = dif(:,K) * (G(K) * invGa(k));
            
            if DEBUG >= 1
                H_sol = y(1:length(PUApprox{1}));
                %PUApprox{1}.sample(H_sol);
                %vol = BlinkVolume(ode,PUApprox{1},t);
                %vol_percent = ((vol-int_vol)/int_vol);
                min_H = min(H_sol);
                fprintf('t = %.5g, delta t = %.3e, min H = %.4g\n',t,h,min_H)
                
                %if min_H<0
                %    warning('H is negative');
                %end
            end
            
            % Predict a solution at t+h.
            tnew = t + h;
            
            if done
                tnew = tfinal;   % Hit end point exactly.
            end
            h = tnew - t;      % Purify h.
            pred = y + sum(dif(:,K),2);
            ynew = pred;
            
            
            % The difference, difkp1, between pred and the final accepted
            % ynew is equal to the backward difference of ynew of order
            % k+1. Initialize to zero for the iteration to compute ynew.
            difkp1 = zeros(neq,1);
            if normcontrol
                normynew = norm(ynew)/length(ynew);
                invwt = 1 / max(max(normy,normynew),threshold);
                minnrm = 100*eps*(normynew * invwt);
            else
                invwt = 1 ./ max(max(abs(y),abs(ynew)),threshold);
                minnrm = 100*eps*norm(ynew .* invwt,inf);
            end
            
            tooslow = false;
            tol_g = 1e-2;
            normres = [];
            for iter = 1:MAXITER    % solve for the proposed ynew value
                
                if use_SNK
                    [rhs,L,U,p,interpnorm] = SNKresidual(tnew,ynew,psi+difkp1,PUApprox,ode,hinvGak,Mtnew,interface_scale,use_parallel,DEBUG);
%                    rhs_interp = NKSresidual(tnew,ynew,hinvGak,PUApprox,ode,interface_scale,use_parallel)...
%						- Masstimes(PUApprox,Mtnew,psi+difkp1);
                else
                    rhs = NKSresidual(tnew,ynew,hinvGak,PUApprox,ode,interface_scale,use_parallel)...
						- Masstimes(PUApprox,Mtnew,psi+difkp1);
                    interpnorm = InterfaceError(PUApprox,rhs_interp)/interface_scale;
                end
                
                normres(iter) = norm(rhs);
                resnorm = normres(iter);
                 
                [lastmsg,lastid] = lastwarn('');
                warning(warnoff);
                
                if iter > 1
                     %tol_g(k) = min(max(abs(normres(k)-linres(k-1))/normres(k-1),tol_g(k-1)^((1+sqrt(5))/2)),1e-2);
                    tol_g(iter) = max(min(tol_g(iter-1),1e-4*(normres(iter)/normres(iter-1))^2),1e-6);
                end
                               
                if use_SNK                    
                    if resnorm==0
                        del = 0;  gmhist = 0;
                    else
                        fun = @(x)SNKjacobian(PUApprox,L,U,p,x,interface_scale,use_parallel,DEBUG);
                        [del,flag,~,~,gmhist] = gmres(fun,-rhs,[],tol_g(iter),100);
                    end
                else
                    [J,L,U,p] = ComputeJacsTime(tnew,ynew,PUApprox,ode,hinvGak,Mtnew,interface_scale);
					fun = @(x)NKSjacobian(PUApprox,J,x,interface_scale,use_parallel);
					prec = @(u)ASPreconditionerTime(PUApprox,L,U,p,u);
                    [del,flag,~,~,gmhist] = gmres(fun,-rhs,[],tol_g(iter),100,prec);
                end
                             
                if DEBUG >= 2
                    fprintf('   GMRES iter = %d, GMRES final = %.3e, resnorm = %.3e, interpnorm = %.3e\n',...
                        length(gmhist)-1,gmhist(end),resnorm,interpnorm)
                    %if length(gmhist) > 50, keyboard, end
                end
                
                warning(warnstat);
                
                % If no new warnings or a muted warning, restore previous lastwarn.
                [msg,msgid] = lastwarn;
                if isempty(msg) || any(strcmp(msgid,warnoffId))
                    lastwarn(lastmsg,lastid);
                end
                
                if normcontrol
                    newnrm = norm(del) * invwt/length(y);
                else
                    newnrm = norm(del .* invwt,inf);
                end
                              
                difkp1 = difkp1 + del;
                ynew = pred + difkp1;
                
                %&& interpnorm<inter_tol && iter>min_iter
                if resnorm<eps || newnrm <= minnrm && interpnorm<inter_tol                   
                    gotynew = true;
                    break;
                elseif iter == 1
                    if havrate
                        errit = newnrm * rate / (1 - rate);
                        if resnorm<eps || errit <= 0.05*rtol && interpnorm<inter_tol     % More stringent when using old rate.
                            gotynew = true;
                            break;
                        end
                    else
                        rate = 0;
                    end
                elseif flag || newnrm > 0.9*oldnrm
                    if iter>min_iter
                        tooslow = true;
                        break;
                    end
                else
                    rate = max(0.9*rate, newnrm / oldnrm);
                    havrate = true;
                    errit = newnrm * rate / (1 - rate);
                    if resnorm<eps || errit <= 0.5*rtol  && interpnorm<inter_tol
                        gotynew = true;
                        break;
                    elseif iter == MAXITER
                        tooslow = true;
                        break;
                    elseif 0.5*rtol < errit*rate^(MAXITER-iter)
                        if iter>min_iter
                            tooslow = true;
                            break;
                        end
                    end
                end
                
                oldnrm = newnrm;
            end                          % end of Newton loop for ynew

            nfevals = nfevals + iter;
            nsolves = nsolves + iter;           
            if tooslow
                nfailed = nfailed + 1;
                if absh <= hmin
                    warning(message('MATLAB:ode15s:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin )));
                    solver_output = odefinalize(solver_name, sol,...
                        outputFcn, outputArgs,...
                        printstats, [nsteps, nfailed, nfevals,...
                        npds, ndecomps, nsolves],...
                        nout, tout, yout,...
                        haveEventFcn, teout, yeout, ieout,...
                        {kvec,dif3d,idxNonNegative});
                    if nargout > 0
                        varargout = solver_output;
                    end
                    return;
                else
                    abshlast = absh;
                    absh = max(0.3 * absh, hmin);
                    h = tdir * absh;
                    done = false;
                    
                    difRU = cumprod((kI - 1 - kJ*(absh/abshlast)) ./ kI) * difU;
                    dif(:,K) = dif(:,K) * difRU(K,K);
                    
                    hinvGak = h * invGa(k);
                    nconhk = 0;
                end
                %        Miter = timeDiff(PUApprox,Mt,dfdy,hinvGak);
                if Mtype == 4
                    Miter = Miter + dMpsidy;
                end
                havrate = false;
            end
        end     % end of while loop for getting ynew
        
        % difkp1 is now the backward difference of ynew of order k+1.
        
        difkp1 = Masstimes(PUApprox,Mtnew,difkp1);
        
        if normcontrol
            err = (norm(difkp1) * invwt) * erconst(k)/length(y);
        else
            err = norm((difkp1) .* invwt,inf) * erconst(k);
        end
        
        if err > rtol                       % Failed step
            nfailed = nfailed + 1;
            if absh <= hmin
                warning(message('MATLAB:ode15s:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin )));
                solver_output = odefinalize(solver_name, sol,...
                    outputFcn, outputArgs,...
                    printstats, [nsteps, nfailed, nfevals,...
                    npds, ndecomps, nsolves],...
                    nout, tout, yout,...
                    haveEventFcn, teout, yeout, ieout,...
                    {kvec,dif3d,idxNonNegative});
                if nargout > 0
                    varargout = solver_output;
                end
                return;
            end
            
            abshlast = absh;
            if nofailed
                nofailed = false;
                hopt = absh * max(0.1, 0.833*(rtol/err)^(1/(k+1))); % 1/1.2
                if k > 1
                    if normcontrol
                        errkm1 = (norm(dif(:,k) + difkp1) * invwt) * erconst(k-1)/length(y);
                    else
                        errkm1 = norm((dif(:,k) + difkp1) .* invwt,inf) * erconst(k-1);
                    end
                    hkm1 = absh * max(0.1, 0.769*(rtol/errkm1)^(1/k)); % 1/1.3
                    if hkm1 > hopt
                        hopt = min(absh,hkm1);      % don't allow step size increase
                        k = k - 1;
                        K = 1:k;
                    end
                end
                absh = max(hmin, hopt);
            else
                absh = max(hmin, 0.5 * absh);
            end
            h = tdir * absh;
            if absh < abshlast
                done = false;
            end
            
            difRU = cumprod((kI - 1 - kJ*(absh/abshlast)) ./ kI) * difU;
            dif(:,K) = dif(:,K) * difRU(K,K);
            
            hinvGak = h * invGa(k);
            nconhk = 0;
            %      Miter = Mt - hinvGak * dfdy;
            if Mtype == 4
                Miter = Miter + dMpsidy;
            end
            havrate = false;
            
        else                                % Successful step
            break;
            
        end
    end % while true
    nsteps = nsteps + 1;
    
    dif(:,k+2) = difkp1 - dif(:,k+1);
    dif(:,k+1) = difkp1;
    for j = k:-1:1
        dif(:,j) = dif(:,j) + dif(:,j+1);
    end
    
    NNreset_dif = false;
    
    %   if haveEventFcn
    %     [te,ye,ie,valt,stop] = odezero(@ntrp15s,eventFcn,eventArgs,valt,...
    %                                    t,y,tnew,ynew,t0,h,dif,k,idxNonNegative);
    %     if ~isempty(te)
    %       if output_sol || (nargout > 2)
    %         teout = [teout, te];
    %         yeout = [yeout, ye];
    %         ieout = [ieout, ie];
    %       end
    %       if stop               % Stop on a terminal event.
    %         % Adjust the interpolation data to [t te(end)].
    %         taux = te(end) - (0:k)*(te(end) - t);
    %         yaux = ntrp15s(taux,t,y,tnew,ynew,h,dif,k,idxNonNegative);
    %         for j=2:k+1
    %           yaux(:,j:k+1) = yaux(:,j-1:k) - yaux(:,j:k+1);
    %         end
    %         dif(:,1:k) = yaux(:,2:k+1);
    %         tnew = te(end);
    %         ynew = ye(:,end);
    %         h = tnew - t;
    %         done = true;
    %       end
    %     end
    %   end
    
    if output_sol
        nout = nout + 1;
        if nout > length(tout)
            tout = [tout, zeros(1,chunk)];  % requires chunk >= refine
            yout = [yout, zeros(neq,chunk)];
            kvec = [kvec, zeros(1,chunk)];
            dif3d = cat(3,dif3d, zeros(neq,maxk+2,chunk));
        end
        tout(nout) = tnew;
        yout(:,nout) = ynew;
        kvec(nout) = k;
        dif3d(:,:,nout) = dif;
    end
    
    if output_ty || haveOutputFcn
        switch outputAt
            case 'SolverSteps'        % computed points, no refinement
                nout_new = 1;
                tout_new = tnew;
                yout_new = ynew;
            case 'RefinedSteps'       % computed points, with refinement
                tref = t + (tnew-t)*S;
                nout_new = refine;
                tout_new = [tref, tnew];
                yout_new = [ntrp15s(tref,[],[],tnew,ynew,h,dif,k,idxNonNegative), ynew];
            case 'RequestedPoints'    % output only at tspan points
                nout_new =  0;
                tout_new = [];
                yout_new = [];
                while next <= ntspan
                    if tdir * (tnew - tspan(next)) < 0
                        if haveEventFcn && stop     % output tstop,ystop
                            nout_new = nout_new + 1;
                            tout_new = [tout_new, tnew];
                            yout_new = [yout_new, ynew];
                        end
                        break;
                    end
                    nout_new = nout_new + 1;
                    tout_new = [tout_new, tspan(next)];
                    if tspan(next) == tnew
                        yout_new = [yout_new, ynew];
                    else
                        yout_new = [yout_new, ntrp15s(tspan(next),[],[],tnew,ynew,h,dif,k,...
                            idxNonNegative)];
                    end
                    next = next + 1;
                end
        end
        
        if nout_new > 0
            if output_ty
                oldnout = nout;
                nout = nout + nout_new;
                if nout > length(tout)
                    tout = [tout, zeros(1,chunk)];  % requires chunk >= refine
                    yout = [yout, zeros(neq,chunk)];
                end
                idx = oldnout+1:nout;
                tout(idx) = tout_new;
                yout(:,idx) = yout_new;
            end
            if haveOutputFcn
                stop = feval(outputFcn,tout_new,yout_new(outputs,:),'',outputArgs{:});
                if stop
                    done = true;
                end
            end
        end
    end
    
    if done
        break
    end
    
    klast = k;
    abshlast = absh;
    nconhk = min(nconhk+1,maxk+2);
    if nconhk >= k + 2
        temp = 1.2*(err/rtol)^(1/(k+1));
        if temp > 0.1
            hopt = absh / temp;
        else
            hopt = 10*absh;
        end
        kopt = k;
        if k > 1
            if normcontrol
                errkm1 = (norm(dif(:,k)) * invwt) * erconst(k-1);
            else
                errkm1 = norm(dif(:,k) .* invwt,inf) * erconst(k-1);
            end
            temp = 1.3*(errkm1/rtol)^(1/k);
            if temp > 0.1
                hkm1 = absh / temp;
            else
                hkm1 = 10*absh;
            end
            if hkm1 > hopt
                hopt = hkm1;
                kopt = k - 1;
            end
        end
        if k < maxk
            if normcontrol
                errkp1 = (norm(dif(:,k+2)) * invwt) * erconst(k+1);
            else
                errkp1 = norm(dif(:,k+2) .* invwt,inf) * erconst(k+1);
            end
            temp = 1.4*(errkp1/rtol)^(1/(k+2));
            if temp > 0.1
                hkp1 = absh / temp;
            else
                hkp1 = 10*absh;
            end
            if hkp1 > hopt
                hopt = hkp1;
                kopt = k + 1;
            end
        end
        if hopt > absh
            absh = hopt;
            if k ~= kopt
                k = kopt;
                K = 1:k;
            end
        end
    end
    
    % Advance the integration one step.
    t = tnew;
    y = ynew;
    if NNreset_dif
        % Used dif for unperturbed solution to select order and interpolate.
        % In perturbing ynew, defined NNidx.  Use now to reset dif to move along
        % constraint.
        dif(NNidx,:) = 0;
    end
    if normcontrol
        normy = normynew;
    end
    Jcurrent = Jconstant;
    switch Mtype
        case {0,1}
            Mcurrent = true;                    % Constant mass matrix I or M.
        case 2
            % M(t) has already been evaluated at tnew in Mtnew.
            Mt = Mtnew;
            Mcurrent = true;
        case {3,4}  % state dependent
            % M(t,y) has not yet been evaluated at the accepted ynew.
            Mcurrent = false;
    end
    
end % while ~done

solver_output = odefinalize(solver_name, sol,...
    outputFcn, outputArgs,...
    printstats, [nsteps, nfailed, nfevals,...
    npds, ndecomps, nsolves],...
    nout, tout, yout,...
    haveEventFcn, teout, yeout, ieout,...
    {kvec,dif3d,idxNonNegative});
if nargout > 0
    varargout = solver_output;
end

end

function solver_output = odefinalize(solver, sol,...
    outfun, outargs,...
    printstats, statvect,...
    nout, tout, yout,...
    haveeventfun, teout, yeout, ieout,...
    interp_data)
%ODEFINALIZE Helper function called by ODE solvers at the end of integration.

%   Jacek Kierzenka
%   Copyright 1984-2005 The MathWorks, Inc.

if ~isempty(outfun)
    feval(outfun,[],[],'done',outargs{:});
end

% Return more stats for implicit solvers: ODE15i, ODE15s, ODE23s, ODE23t, ODE23tb
fullstats = (length(statvect) > 3);  % faster than 'switch' or 'ismember'

stats = struct('nsteps',statvect(1),'nfailed',statvect(2),'nfevals',statvect(3));
if fullstats
    stats.npds     = statvect(4);
    stats.ndecomps = statvect(5);
    stats.nsolves  = statvect(6);
else
    statvect(4:6) = 0;   % Backwards compatibility
end

if printstats
    fprintf(getString(message('MATLAB:odefinalize:LogSuccessfulSteps', sprintf('%g',stats.nsteps))));
    fprintf(getString(message('MATLAB:odefinalize:LogFailedAttempts', sprintf('%g',stats.nfailed))));
    fprintf(getString(message('MATLAB:odefinalize:LogFunctionEvaluations', sprintf('%g',stats.nfevals))));
    if fullstats
        fprintf(getString(message('MATLAB:odefinalize:LogPartialDerivatives', sprintf('%g',stats.npds))));
        fprintf(getString(message('MATLAB:odefinalize:LogLUDecompositions', sprintf('%g',stats.ndecomps))));
        fprintf(getString(message('MATLAB:odefinalize:LogSolutionsOfLinearSystems', sprintf('%g',stats.nsolves))));
    end
end

solver_output = {};

if (nout > 0) % produce output
    if isempty(sol) % output [t,y,...]
        solver_output{1} = tout(1:nout).';
        solver_output{2} = yout(:,1:nout).';
        if haveeventfun
            solver_output{3} = teout.';
            solver_output{4} = yeout.';
            solver_output{5} = ieout.';
        end
        solver_output{end+1} = statvect(:);  % Column vector
    else % output sol
        % Add remaining fields
        sol.x = tout(1:nout);
        sol.y = yout(:,1:nout);
        if haveeventfun
            sol.xe = teout;
            sol.ye = yeout;
            sol.ie = ieout;
        end
        sol.stats = stats;
        [kvec,dif3d,idxNonNegative] = deal(interp_data{:});
        sol.idata.kvec = kvec(1:nout);
        maxkvec = max(sol.idata.kvec);
        sol.idata.dif3d = dif3d(:,1:maxkvec+2,1:nout);
        sol.idata.idxNonNegative = idxNonNegative;
        solver_output{1} = sol;
    end
end

end

function [Jconstant,Jfcn,Jargs,Joptions] = ...
    odejacobian(fcnHandlesUsed,ode,t0,y0,options,extras)
%ODEJACOBIAN  Helper function for the Jacobian function in ODE solvers
%    ODEJACOBIAN determines whether the Jacobian is constant and if so,
%    returns its value as Jfcn. If an analytical Jacobian is available from
%    a function, ODEJACOBIAN initializes Jfcn and creates a cell array of
%    additional input arguments. For numerical Jacobian, ODEJACOBIAN tries to
%    extract JPattern and sets JOPTIONS for use with ODENUMJAC.
%
%   See also ODE15S, ODE23S, ODE23T, ODE23TB, ODENUMJAC.

%   Jacek Kierzenka
%   Copyright 1984-2009 The MathWorks, Inc.

Jconstant = strcmp(odeget(options,'JConstant','off','fast'),'on');
Jfcn = [];
Jargs = {};
Joptions = [];

Janalytic = false;

if fcnHandlesUsed
    Jfcn = odeget(options,'Jacobian',[],'fast');
    if ~isempty(Jfcn)
        if isnumeric(Jfcn)
            Jconstant = true;
        else
            Janalytic = true;
            Jargs = extras;
        end
    end
else  % ode-file used
    joption = odeget(options,'Jacobian','off','fast');
    switch lower(joption)
        case 'on'    % ode(t,y,'jacobian',p1,p2...)
            Janalytic = true;
            Jfcn = ode;
            Jargs = [{'jacobian'} extras];
        case 'off'   % use odenumjac
        otherwise
            error(message('MATLAB:odejacobian:InvalidJOption', joption));
    end
end

if ~Janalytic   % odenumjac will be used
    Joptions.diffvar  = 2;       % df(t,y)/dy
    Joptions.vectvars = [];
    vectorized = strcmp(odeget(options,'Vectorized','off','fast'),'on');
    if vectorized
        Joptions.vectvars = 2;     % f(t,[y1,y2]) = [f(t,y1), f(t,y2)]
    end
    
    atol = odeget(options,'AbsTol',1e-6,'fast');
    Joptions.thresh = zeros(size(y0))+ atol(:);
    Joptions.fac  = [];
    
    if fcnHandlesUsed
        jpattern = odeget(options,'JPattern',[],'fast');
    else  % ode-file used
        jp_option = odeget(options,'JPattern','off','fast');
        switch lower(jp_option)
            case 'on'
                jpattern = feval(ode,[],[],'jpattern',extras{:});
            case 'off'  % no pattern provided
                jpattern = [];
            otherwise
                error(message('MATLAB:odejacobian:InvalidJpOption', jp_option));
        end
    end
    if ~isempty(jpattern)
        Joptions.pattern = jpattern;
        Joptions.g = colgroup(jpattern);
    end
end

end

function [neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, args, odeFcn, ...
    options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, ...
    dataType ] =   ...
    odearguments(FcnHandlesUsed, solver, ode, tspan, y0, options, extras)
%ODEARGUMENTS  Helper function that processes arguments for all ODE solvers.
%
%   See also ODE113, ODE15I, ODE15S, ODE23, ODE23S, ODE23T, ODE23TB, ODE45.

%   Mike Karr, Jacek Kierzenka
%   Copyright 1984-2015 The MathWorks, Inc.

if strcmp(solver,'ode15i')
    FcnHandlesUsed = true;   % no MATLAB v. 5 legacy for ODE15I
end

if FcnHandlesUsed  % function handles used
    if isempty(tspan) || isempty(y0)
        error(message('MATLAB:odearguments:TspanOrY0NotSupplied', solver));
    end
    if length(tspan) < 2
        error(message('MATLAB:odearguments:SizeTspan', solver));
    end
    htspan = abs(tspan(2) - tspan(1));
    tspan = tspan(:);
    ntspan = length(tspan);
    t0 = tspan(1);
    next = 2;       % next entry in tspan
    tfinal = tspan(end);
    args = extras;                 % use f(t,y,p1,p2...)
    
else  % ode-file used   (ignored when solver == ODE15I)
    % Get default tspan and y0 from the function if none are specified.
    if isempty(tspan) || isempty(y0)
        if exist(ode)==2 && ( nargout(ode)<3 && nargout(ode)~=-1 )
            error(message('MATLAB:odearguments:NoDefaultParams', funstring( ode ), solver, funstring( ode )));
        end
        [def_tspan,def_y0,def_options] = feval(ode,[],[],'init',extras{:});
        if isempty(tspan)
            tspan = def_tspan;
        end
        if isempty(y0)
            y0 = def_y0;
        end
        options = odeset(def_options,options);
    end
    tspan = tspan(:);
    ntspan = length(tspan);
    if ntspan == 1    % Integrate from 0 to tspan
        t0 = 0;
        next = 1;       % Next entry in tspan.
    else
        t0 = tspan(1);
        next = 2;       % next entry in tspan
    end
    htspan = abs(tspan(next) - t0);
    tfinal = tspan(end);
    
    % The input arguments of f determine the args to use to evaluate f.
    %   if (exist(ode)==2)
    %     if (nargin(ode) == 2)
    %       args = {};                   % f(t,y)
    %     else
    %       args = [{''} extras];        % f(t,y,'',p1,p2...)
    %     end
    %   else  % MEX-files, etc.
    %     try
    %       args = [{''} extras];        % try f(t,y,'',p1,p2...)
    %       feval(ode,tspan(1),y0(:),args{:});
    %     catch
    %       args = {};                   % use f(t,y) only
    %     end
    %   end
    args = {};
end

y0 = y0(:);
neq = length(y0);

% Test that tspan is internally consistent.
if any(isnan(tspan))
    error(message('MATLAB:odearguments:TspanNaNValues'));
end
if t0 == tfinal
    error(message('MATLAB:odearguments:TspanEndpointsNotDistinct'));
end
tdir = sign(tfinal - t0);
if any( tdir*diff(tspan) <= 0 )
    error(message('MATLAB:odearguments:TspanNotMonotonic'));
end

%f0 = feval(ode,t0,y0,args{:});   % ODE15I sets args{1} to yp0.
f0 = y0;

[m,n] = size(f0);
if n > 1
    error(message('MATLAB:odearguments:FoMustReturnCol', funstring( ode )));
elseif m ~= neq
    error(message('MATLAB:odearguments:SizeIC', funstring( ode ), m, neq, funstring( ode )));
end

% Determine the dominant data type
classT0 = class(t0);
classY0 = class(y0);
classF0 = class(f0);
if strcmp(solver,'ode15i')
    classYP0 = class(args{1});  % ODE15I sets args{1} to yp0.
    dataType = superiorfloat(t0,y0,args{1},f0);
    
    if ~( strcmp(classT0,dataType) && strcmp(classY0,dataType) && ...
            strcmp(classF0,dataType) && strcmp(classYP0,dataType))
        input1 = '''t0'', ''y0'', ''yp0''';
        input2 = '''f(t0,y0,yp0)''';
        warning(message('MATLAB:odearguments:InconsistentDataType',input1,input2,solver));
    end
else
    dataType = superiorfloat(t0,y0,f0);
    
    if ~( strcmp(classT0,dataType) && strcmp(classY0,dataType) && ...
            strcmp(classF0,dataType))
        input1 = '''t0'', ''y0''';
        input2 = '''f(t0,y0)''';
        warning(message('MATLAB:odearguments:InconsistentDataType',input1,input2,solver));
    end
end

% Get the error control options, and set defaults.
rtol = odeget(options,'RelTol',1e-3,'fast');
if (length(rtol) ~= 1) || (rtol <= 0)
    error(message('MATLAB:odearguments:RelTolNotPosScalar'));
end
if rtol < 100 * eps(dataType)
    rtol = 100 * eps(dataType);
    warning(message('MATLAB:odearguments:RelTolIncrease', sprintf( '%g', rtol )))
end
atol = odeget(options,'AbsTol',1e-6,'fast');
if any(atol <= 0)
    error(message('MATLAB:odearguments:AbsTolNotPos'));
end
normcontrol = strcmp(odeget(options,'NormControl','off','fast'),'on');
if normcontrol
    if length(atol) ~= 1
        error(message('MATLAB:odearguments:NonScalarAbsTol'));
    end
    normy = norm(y0);
else
    if (length(atol) ~= 1) && (length(atol) ~= neq)
        error(message('MATLAB:odearguments:SizeAbsTol', funstring( ode ), neq));
    end
    atol = atol(:);
    normy = [];
end
threshold = atol / rtol;

% By default, hmax is 1/10 of the interval.
hmax = min(abs(tfinal-t0), abs(odeget(options,'MaxStep',0.1*(tfinal-t0),'fast')));
if hmax <= 0
    error(message('MATLAB:odearguments:MaxStepLEzero'));
end
htry = abs(odeget(options,'InitialStep',[],'fast'));
if ~isempty(htry) && (htry <= 0)
    error(message('MATLAB:odearguments:InitialStepLEzero'));
end

odeFcn = ode;

end

function [haveeventfun,eventFcn,eventArgs,eventValue,teout,yeout,ieout] =...
    odeevents(FcnHandlesUsed,ode,t0,y0,options,extras)
%ODEEVENTS  Helper function for the events function in ODE solvers
%    ODEEVENTS initializes eventFcn to the events function, and creates a
%    cell-array of its extra input arguments. ODEEVENTS evaluates the events
%    function at(t0,y0).
%
%   See also ODE113, ODE15S, ODE23, ODE23S, ODE23T, ODE23TB, ODE45.

%   Jacek Kierzenka
%   Copyright 1984-2010 The MathWorks, Inc.

haveeventfun = 0;   % no Events function
eventArgs = [];
eventValue = [];
teout = [];
yeout = [];
ieout = [];

eventFcn = odeget(options,'Events',[],'fast');
if isempty(eventFcn)
    return
end

if FcnHandlesUsed     % function handles used
    haveeventfun = 1;   % there is an Events function
    eventArgs = extras;
    eventValue = feval(eventFcn,t0,y0,eventArgs{:});
    
else   % ode-file used
    switch lower(eventFcn)
        case 'on'
            haveeventfun = 1;   % there is an Events function
            eventFcn = ode;            % call ode(t,y,'events',p1,p2...)
            eventArgs = [{'events'}, extras];
            eventValue = feval(eventFcn,t0,y0,eventArgs{:});
        case 'off'
        otherwise
            error(message('MATLAB:odeevents:MustSetOnOrOff'))
    end
    
end

end