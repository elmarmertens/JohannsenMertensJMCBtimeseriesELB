function results=ols(y, X, Tscaleflag)
% PURPOSE: least-squares regression 
%---------------------------------------------------
% USAGE: results = ols(y,X, Tscaleflag)
% where: y = dependent variable vector    (nobs x 1)
%        X = independent variables matriX (nobs x nvar)
%        Tscaleflag = flag. If set to one uses MLE of sigma_e (like cov(x,1))
%---------------------------------------------------
% RETURNS: a structure
%        results.meth  = 'ols'
%        results.Tscale= Tscaleflag
%        results.beta  = bhat     (nvar x 1)
%        results.Vbeta = Vhat-matrix of bhat     (nvar x nvar)
%        results.tstat = t-stats  (nvar x 1)
%        results.yhat  = yhat     (nobs x 1)
%        results.resid = residuals (nobs x 1)
%        results.ssr   = e'*e;
%        results.sige  = e'*e/(n-k)   scalar
%        results.rsqr  = rsquared     scalar
%        results.rbar  = rbar-squared scalar
%        results.dw    = Durbin-Watson Statistic
%        results.nobs  = nobs
%        results.nvar  = nvars
%        results.y     = y data vector (nobs x 1)
%        results.X     = X data vector (nobs x nvar)
%        results.bint  = (nvar x 2) vector with 95% confidence intervals on beta
%---------------------------------------------------
% SEE ALSO: prt(results), plt(results)
%---------------------------------------------------

%   Coded by  Elmar Mertens, em@elmarmertens.com

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com
%
% Barry Dillon (CICG Equity)
% added the 95% confidence intervals on bhat
% 
% amendments by:
% Elmar Mertens
% Graduate Student and TA at the University of Basel, 
% WWZ/Department of Finance
% email: elmar.mertens@unibas.ch
% added SSR, X data field, Vbeta and Tscaleflag

if (nargin < 2) 
   error('Wrong # of arguments'); 
elseif (nargin == 2)
   Tscaleflag = 0;
elseif (nargin >3)
   error('Wrong # of arguments'); 
end

[nobs nvar] = size(X); [nobs2 junk] = size(y);

if (nobs ~= nobs2); error('X and y must have same # obs in ols');  end;

results.meth = 'ols';
results.Tscale= Tscaleflag;
results.y = y;
results.X = X;
results.nobs = nobs;
results.nvar = nvar;

[q r] = qr(X,0);
XpXi = (r'*r)\eye(nvar);


results.beta = XpXi*(X'*y);
results.yhat = X*results.beta;
results.resid = y - results.yhat;
results.ssr = results.resid'*results.resid;
if Tscaleflag
   results.sige = results.ssr/nobs;
else
   results.sige = results.ssr/(nobs-nvar);
end
results.Vbeta = (results.sige)*(XpXi);
sigb=sqrt(diag(results.Vbeta));
tcrit=-tdis_inv(.025,nobs);
results.bint=[results.beta-tcrit.*sigb, results.beta+tcrit.*sigb];
results.tstat = results.beta./(sigb);
if all(y == 1)
   [results.rsqr, results.rbar, results.dw] = deal(NaN);
else
   ym = y - mean(y);
   rsqr1 = results.ssr;
   rsqr2 = ym'*ym;
   results.rsqr = 1.0 - rsqr1/rsqr2; % r-squared
   rsqr1 = rsqr1/(nobs-nvar);
   rsqr2 = rsqr2/(nobs-1.0);
   if rsqr2 ~= 0
      results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
   else
      results.rbar = results.rsqr;
   end;
   ediff = results.resid(2:nobs) - results.resid(1:nobs-1);
   results.dw = (ediff'*ediff)/results.ssr; % durbin-watson
end 

results.llf = - nobs / 2 * (1 + log(2 * pi) + log(results.ssr / nobs)); %for compatibility with EViews I do not check Tscaleflag here
