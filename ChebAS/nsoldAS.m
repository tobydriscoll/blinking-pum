function [sol, it_hist, ierr, x_hist] = nsold(x,evalF,tol,num_sols,PUApprox,parms)
% NSOLDAS  Newton-Armijo nonlinear solver
%
% Factor Jacobians with Gaussian Elimination
%
% Hybrid of Newton, Shamanskii, Chord
%
% C. T. Kelley, April 1, 2003.
%
% This code comes with no guarantee or warranty of any kind.
%
% function [sol, it_hist, ierr, x_hist] = nsold(x,f,tol,parms)
%
% inputs:
%        initial iterate = x
%        function = f
%        tol = [atol, rtol] relative/absolute
%                           error tolerances
%        parms = [maxit, isham, rsham, jdiff, nl, nu]
%        maxit = maxmium number of iterations
%                default = 40
%        isham, rsham: The Jacobian matrix is
%                computed and factored after isham
%                updates of x or whenever the ratio
%                of successive l2 norms of the
%                 nonlinear residual exceeds rsham.
%
%            isham = -1, rsham = .5 is the default
%            isham =  1, rsham = 0 is Newton's method,
%            isham = -1, rsham = 1 is the chord method,
%            isham =  m, rsham = 1 is the Shamanskii method with
%                        m steps per Jacobian evaluation
%
%                       The Jacobian is computed and factored
%                       whenever the stepsize
%                       is reduced in the line search.
%
%       jdiff = 1: compute Jacobians with forward differences
%       jdiff = 0: a call to f will provide analytic Jacobians
%                         using the syntax [function,jacobian] = f(x)
%                 defaults = [40, 1000, .5, 1]
%
%       nl, nu: lower and upper bandwidths of a banded Jacobian.
%               If you include nl and nu in the parameter list,
%               the Jacobian will be evaluated with a banded differencing
%               scheme and stored as a sparse matrix.
%
%
%
% output:
%    sol = solution
%    it_hist = array of iteration history, useful for tables and plots
%                The two columns are the residual norm and
%                number of step size reductions done in the line search.
%
%        ierr = 0 upon successful termination
%        ierr = 1 if after maxit iterations
%             the termination criterion is not satsified
%        ierr = 2 failure in the line search. The iteration
%             is terminated if too many steplength reductions
%             are taken.
%
%    x_hist = matrix of the entire interation history.
%             The columns are the nonlinear iterates. This
%             is useful for making movies, for example, but
%             can consume way too much storage. This is an
%             OPTIONAL argument. Storage is only allocated
%             if x_hist is in the output argument list.
%
%
% internal parameter:
%       debug = turns on/off iteration statistics display as
%               the iteration progresses
%
% Here is an example. The example computes pi as a root of sin(x)
% with Newton's method, forward difference derivatives,
% and plots the iteration history. Note that x_hist is not in
% the output list.
%
%
%  x = 3; tol = [1.d-6, 1.d-6]; params = [40, 1, 0];
%  [result, errs, ierr] = nsold(x, 'sin', tol, params);
%  result
%  semilogy(errs)
%
%
% Set the debug parameter, 1 turns display on, otherwise off.
%
debug = 0;
%
% Initialize it_hist, ierr, and set the iteration parameters.
%
ierr = 0;
maxarm = 20;
maxit = 40;
isham = -1;
rsham = .5;
jdiff = 1;
iband = 0;


if nargin >= 6 & length(parms) ~= 0
    maxit = parms(1); isham = parms(2); rsham = parms(3);
    if length(parms) >= 4
        jdiff = parms(4);
    end
    if length(parms) >= 6
        nl = parms(5); nu = parms(6);
        iband = 1;
    end
end

rtol = tol(2); atol = tol(1);
it_hist = [];
n = length(x);
if nargout == 4, x_hist = x; end
fnrm = 1;
itc = 0;
%
% evaluate f at the initial iterate
% compute the stop tolerance
%
%f0 = feval(f,x);

%precompute interpolation matricies
setInterpMatrices(PUApprox);

f0 = ParPreconditionedNewtonForward(PUApprox,x,evalF,num_sols);

fnrm = norm(f0);
it_hist = [fnrm,0];
fnrmo = 1;
itsham = isham;
stop_tol = atol+rtol*fnrm;
%
% main iteration loop
%
while(fnrm > stop_tol & itc < maxit)
%
% keep track of the ratio (rat = fnrm/frnmo)
% of successive residual norms and 
% the iteration counter (itc)
%
    rat = fnrm/fnrmo;
    outstat(itc+1, :) = [itc fnrm rat];
    fnrmo = fnrm; 
    itc = itc+1;
%
% evaluate and factor the Jacobian
% on the first iteration, every isham iterates, or
% if the ratio of successive residual norm is too large
%
    if(itc == 1 | rat > rsham | itsham == 0 | armflag == 1)
        itsham = isham;
    jac_age = -1;
            %Replace this with our new cool solver
            %[fv,jac] = feval(f,x);
            [ fv,jac ] = ParPreconditionedNewtonForward(PUApprox,x,evalF,num_sols);
    end
    itsham = itsham-1;
%
% compute the Newton direction 
%

    %Replace this with krylov method
    %tmp = -l\f0;
    %direction = u\tmp;
    
    [direction] = gmres(@(x)JacobianFoward(PUApprox,jac,x),-f0,[],1e-4,100,[],[],zeros(size(x)));
%
% Add one to the age of the Jacobian after the factors have been
% used in a solve. A fresh Jacobian has an age of -1 at birth.
%
    jac_age = jac_age+1;
    xold = x; fold = f0;
    [step,iarm,x,f0,armflag] = armijo(direction,x,f0,evalF,maxarm,PUApprox,num_sols);
%
% If the line search fails and the Jacobian is old, update it.
% If the Jacobian is fresh; you're dead.
%
    if armflag == 1  
       if jac_age > 0
          sol = xold;
          x = xold; f0 = fold;    
          disp('Armijo failure; recompute Jacobian.');
       else
          disp('Complete Armijo failure.');
          sol = xold;
          ierr = 2;
          return
       end
    end
    
    
    fnrm = norm(f0);
    it_hist = [it_hist',[fnrm,iarm]']';
    if nargout == 4, x_hist = [x_hist,x]; end
    rat = fnrm/fnrmo;
    if debug == 1, disp([itc fnrm rat]); end
    outstat(itc+1, :) = [itc fnrm rat];
% end while
end
sol = x;
if debug == 1, disp(outstat); end
%
% on failure, set the error flag
%
if fnrm > stop_tol, ierr = 1; end
%
%
%

function [step,iarm,xp,fp,armflag] = armijo(direction,x,f0,evalF,maxarm,PUApprox,num_sols)
iarm = 0;
sigma1 = .5;
alpha = 1.d-4;
armflag = 0;
xp = x; fp = f0; 
%
    xold = x;
    lambda = 1; lamm = 1; lamc = lambda; iarm = 0;
    step = lambda*direction;
    xt = x + step;
    
    %ft = feval(f,xt);
    ft = ParPreconditionedNewtonForward(PUApprox,xt,evalF,num_sols);
    
    nft = norm(ft); nf0 = norm(f0); ff0 = nf0*nf0; ffc = nft*nft; ffm = nft*nft;
    while nft >= (1 - alpha*lambda) * nf0;
%
%   Apply the three point parabolic model.
%
        if iarm == 0
            lambda = sigma1*lambda;
        else
            lambda = parab3p(lamc, lamm, ff0, ffc, ffm);
        end
%
% Update x; keep the books on lambda.
%
        step = lambda*direction;
        xt = x + step;
        lamm = lamc;
        lamc = lambda;
%
% Keep the books on the function norms.
%
        %ft = feval(f,xt);
        ft = ParPreconditionedNewtonForward(PUApprox,ft,evalF,num_sols);
        nft = norm(ft);
        ffm = ffc;
        ffc = nft*nft;
        iarm = iarm+1;
        if iarm > maxarm
            disp(' Armijo failure, too many reductions ');
            armflag = 1;
            sol = xold;
            return;
        end
    end
    xp = xt; fp = ft;
%
%   end of line search
%

%
function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
% Apply three-point safeguarded parabolic model for a line search.
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
%
% input:
%       lambdac = current steplength
%       lambdam = previous steplength
%       ff0 = value of \| F(x_c) \|^2
%       ffc = value of \| F(x_c + \lambdac d) \|^2
%       ffm = value of \| F(x_c + \lambdam d) \|^2
%
% output:
%       lambdap = new value of lambda given parabolic model
%
% internal parameters:
%       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
%

%
% Set internal parameters.
%
sigma0 = .1; sigma1 = .5;
%
% Compute coefficients of interpolation polynomial.
%
% p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1
%
% d1 = (lambdac - lambdam)*lambdac*lambdam < 0
%      so, if c2 > 0 we have negative curvature and default to
%      lambdap = sigam1 * lambda.
%
c2 = lambdam*(ffc-ff0)-lambdac*(ffm-ff0);
if c2 >= 0
    lambdap = sigma1*lambdac; return
end
c1 = lambdac*lambdac*(ffm-ff0)-lambdam*lambdam*(ffc-ff0);
lambdap = -c1*.5/c2;
if lambdap < sigma0*lambdac, lambdap = sigma0*lambdac; end
if lambdap > sigma1*lambdac, lambdap = sigma1*lambdac; end