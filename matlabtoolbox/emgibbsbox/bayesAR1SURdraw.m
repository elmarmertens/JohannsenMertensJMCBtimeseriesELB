function [rhodraw, resid, Erho, Vrho] = bayesAR1SURdraw(y, ylag, Vresid, rho0, rhoV0i, Ndraws, rndStream)
% BAYESREGRESSIONSLOPESGIBBSDRAW performs Gibbs steps for linear regression model 
%   this function draws slopes conditional on precision (inverse residual variance)
%   Usage: [bdraw, resid, b, V] = bayesRegressionSlopesGibbsDraw(y, X, h, b0, V0, Ndraws, rndStream)
%
% See also bayesRegressionGibbsDraw 

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 17-Jan-2009 16:30:05 $
% $Revision : 1.00 $
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : bayesRegressionGibbsDraw.m

if nargin < 6 || isempty(Ndraws)
   Ndraws = 1;
end
if nargin < 7 || isempty(rndStream)
   rndStream = RandStream.getDefaultStream;
end


[T, Ny] = size(y);
resid   = NaN(T,Ny);
% rhodraw = NaN(Ny,1);
Iy      = eye(Ny);

H      = Iy / Vresid;
ytilde = y * H;
XX     = ylag' * ylag;

rhoVi = rhoV0i + XX .* H;
Vrho  = Iy / rhoVi;
Erho  = Vrho * (rhoV0i * rho0 + sum(ylag .* ytilde)');

cholVrho = chol(Vrho)';
rhodraw = bsxfun(@plus, Erho, cholVrho * randn(rndStream, Ny, Ndraws)); 

if nargout > 1 && (Ndraws == 1)
   resid       = y - bsxfun(@times, ylag, rhodraw');
end
   
