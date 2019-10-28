function [bdraw, resid, b, V] = bayesRegressionSlopesGibbsDrawDiffuse(y, X, h, Ndraws, rndStream)
% BAYESREGRESSIONSLOPESGIBBSDRAW performs Gibbs steps for linear regression model 
%   this function draws slopes conditional on precision (inverse residual variance)
%   Usage: [bdraw, resid, b, V] = bayesRegressionSlopesGibbsDraw(y, X, h, Ndraws, rndStream)
%
% See also bayesRegressionGibbsDraw 

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 17-Jan-2009 16:30:05 $
% $Revision : 1.00 $
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : bayesRegressionGibbsDraw.m

if nargin < 4 || isempty(Ndraws)
   Ndraws = 1;
end
% if nargin < 5 || isempty(rndStream)
%    rndStream = RandStream.getDefaultStream;
% end

%% Some prelim transformations
[T, K]      = size(X); %#ok<ASGLU>
I           = eye(K);

if nargin < 4 || isempty(Ndraws)
   Ndraws = 1;
end
XX    = X' * X;
Xy    = X' * y;


%% draw Regression-slope

% Note: 
% The choleski of V will be needed anyway for the random number generation, 
% hence it is (slightly) more efficient to use the choleski also for computing the inverse

Vi    = XX * h; 
cV    = chol(Vi) \ I; % notice: Matlab's choleski delivers UPPER triangular matrix
V     = cV * cV';

b     = V * (Xy * h);

bdraw = bsxfun(@plus, b, cV * randn(rndStream, size(b,1), Ndraws)); 

if nargout > 1
    resid       = bsxfun(@minus, y, X * bdraw);
end
   
