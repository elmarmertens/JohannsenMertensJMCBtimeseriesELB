function results=nwest(y,X,nlag)
% PURPOSE: computes Newey-West adjusted heteroscedastic-serial
%          consistent Least-squares Regression, including constant
%---------------------------------------------------
% USAGE: results = nwest(y,X,nlag)
% where: y = dependent variable vector (nobs x 1)
%        X = independent variables matrix (nobs x nvar)
%     nlag = lag length to use
%---------------------------------------------------
% RETURNS: a structure
%        results.meth  = 'newlyw'
%        results.beta  = bhat
%        results.tstat = t-stats
%        results.yhat  = yhat
%        results.resid = residuals
%        results.sige  = e'*e/(n-k)
%        results.rsqr  = rsquared
%        results.rbar  = rbar-squared
%        results.dw    = Durbin-Watson Statistic
%        results.nobs  = nobs
%        results.nvar  = nvars
%        results.y     = y data vector
% --------------------------------------------------
% SEE ALSO: nwest_d, prt(results), plt(results)
%---------------------------------------------------
% References:  Gallant, R. (1987),
%  "Nonlinear Statistical Models," pp.137-139.
%---------------------------------------------------


% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% % jlesage@spatial-econometrics.com

% Amendments made by Elmar Mertens, em@elmarmertens.com
% - nwlags get stored, optimal lag-length for nlag by default

error(nargchk(2,3,nargin))
[nobs nvar] = size(X);

if nargin < 3 || isempty(nlag)
    nlag = floor(4 * (nobs / 100)^(2/9));
end


results.meth    = 'nwest';
results.y       = y;
results.X       = X;
results.nobs    = nobs;
results.nvar    = nvar;
results.nlag    = nlag;

% xpxi            = inv(X'*X);
[q, r] = qr(X,0);
xpxi   = (r'*r)\eye(nvar);

results.beta    = xpxi*(X'*y);
results.yhat    = X*results.beta;
results.resid   = y - results.yhat;
sigu = results.resid'*results.resid;
results.sige    = sigu/(nobs-nvar);

% perform Newey-West correction

        %V=xpxi*G*xpxi;
        orth = results.resid(:, ones(nvar, 1)) .* X;
        %V = xpxi * haccme(orth', nlag) * xpxi * nobs^2 / (nobs - nvar);
        V = xpxi * haccme(orth', nlag) * xpxi * nobs;
        nwerr= sqrt(diag(V));

results.Vbeta = V;
results.tstat = results.beta./nwerr; % Newey-West t-statistics
ym = y - ones(nobs,1)*mean(y);
rsqr1 = sigu;
rsqr2 = ym'*ym;
results.rsqr = 1.0 - rsqr1/rsqr2; % r-squared
rsqr1 = rsqr1/(nobs-nvar);
rsqr2 = rsqr2/(nobs-1.0);
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
ediff = results.resid(2:nobs) - results.resid(1:nobs-1);
results.dw = diag((ediff'*ediff)./(sigu))'; % durbin-watson

