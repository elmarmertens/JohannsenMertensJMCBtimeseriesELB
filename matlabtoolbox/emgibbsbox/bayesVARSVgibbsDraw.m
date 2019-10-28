function [aDraws, resid, kompanion, maxlambda, eigenvalues, eigenvectors, a, aSigma]  = bayesVARSVgibbsDraw(Y, p, constantFlag, iSigmaResid, a0, aSigma0, rndStream, rejectionTrials, one)
% bayesVARgibbsDraw performs Gibbs step for vector linear regression model with known variance
%   This Gibbs step draws slopes conditional on the VCV matrix of the residuals
% Usage:
%
%   [aDraws, resid, kompanion, maxlambda, a, aSigma]  = ...
%       bayesVARSVgibbsDraw(Y, p, constantFlag, iSigmaResid, a0, aSigma0, ...
%       rndStream, rejectionTrials, one)
%
%   Y = X * A + resid, where A = reshape(aDraw, Nx, Ny)
%
%   and 
%
%     varcoeff  = reshape(aDraws(:, trial), Nx, Ny)';
%     if constantFlag
%        varcoeff = varcoeff(:,1:end-1);
%     end
%     kompanion         = [varcoeff; eye(Ny*(p-1), Ny*p)];
%       
%
%   [aDraws, resid, kompanion, maxlambda, a, aSigma] = bayesVARgibbsDraw(Y, p, ...
%     constantFlag, iSigmaResid, a0, aSigma0, rndStream, rejectionTrials)
%  
%  or
% 
%   [aDraws, ...] = bayesVARgibbsDraw(Y, X, ...);
% 
% See also bayesVectorRegressionGibbsDraw2Steps, bayesVectorRegressionGibbsDraw, lag4VAR

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 01-Sep-2010 11:35:11 $
% $Revision : 1.00 $
% DEVELOPED : 7.8.0.347 (R2009a)
% FILENAME  : bayesVARgibbsDraw.m

if isscalar(p)
   [X, Y]   = lag4VAR(Y, p, constantFlag);
else
   X = p;
   p = NaN;
end

if nargin < 7 || isempty(rndStream)
   rndStream = getDefaultStream;
end

if nargin < 8 || isempty(rejectionTrials)
   rejectUnstable = false;
   Ndraw          = 1;
else
   rejectUnstable = true;
   Ndraw          = rejectionTrials;
end

if nargin < 9 || isempty(one)
    one = 1;
end

if nargout > 4
    reportEigenvalues = true;
else
    reportEigenvalues = false;
end

[T, Ny]  = size(Y);
Nx       = size(X, 2);
Na       = length(a0);
Ia       = eye(Na);
% Iy       = eye(Ny);

%% posterior for coefficients

% posterior variance
inv_aSigma0    = Ia / aSigma0;
inv_aSigma     = inv_aSigma0;
for t = 1 : T
    inv_aSigma = inv_aSigma + kron(iSigmaResid(:,:,t), X(t,:)' * X(t,:));
end
chol_aSigma    = chol(inv_aSigma) \ Ia;            % notice: Matlab's choleski delivers UPPER triangular matrix
aSigma         = chol_aSigma * chol_aSigma';       % checkdiff(aSigma, inv(inv_aSigma));

% posterior mean
Ytilde = zeros(size(Y));
for t = 1 : T
    Ytilde(t,:) = Y(t,:) * iSigmaResid(:,:,t);
end
tmp      = X' * Ytilde;
a        = aSigma * (inv_aSigma0 * a0 + tmp(:));

% draw from posterior
aDraws   = bsxfun(@plus, a, chol_aSigma * randn(rndStream, Na, Ndraw)); % chol_aSigma is the UPPER triangular factorization of aSigma, but this is OK for drawing RV

if rejectUnstable
   
   % need to know p, in order to construct companion
   if ~exist('p', 'var')
      NxOhneKonstante = Nx - double(constantFlag);
      p = NxOhneKonstante / Ny;
      % check
      if mod(NxOhneKonstante, Ny) ~= 0
         error('NxOhneKonstant should be divisble by Ny, but we have NxOhneKonstant=%d and Ny=%d', NxOhneKonstant, Ny)
      end
   end
   
   stabledraw = false;
   trial = 0;
   while ~stabledraw && trial < rejectionTrials
      trial         = trial + 1;
      
      varcoeff          = reshape(aDraws(:, trial), Nx, Ny);
      if constantFlag
         varcoeff = varcoeff(1:end-1,:);
      end
      kompanion         = [varcoeff'; eye(Ny*(p-1), Ny*p)];
      if reportEigenvalues
          [eigenvectors, Lambda] = eig(kompanion);
          eigenvalues = diag(Lambda);
          maxlambda         = max(abs(eigenvalues));
      else
          maxlambda         = max(abs(eig(kompanion)));
      end
      
      stabledraw        = maxlambda < one;
      
   end
   if ~stabledraw
      error('No stable draw could be found after %d trials', trial);
   end
   
   aDraws = aDraws(:,trial);
   
elseif nargout > 2
   varcoeff          = reshape(aDraws, Nx, Ny);
   if constantFlag
      varcoeff = varcoeff(1:end-1,:);
   end
   kompanion         = [varcoeff'; eye(Ny*(p-1), Ny*p)];
   if nargout > 3
       if reportEigenvalues
           [eigenvectors, Lambda] = eig(kompanion);
           eigenvalues = diag(Lambda);
           maxlambda         = max(abs(eigenvalues));
       else
           maxlambda         = max(abs(eig(kompanion)));
       end
   end
end

if nargout > 1
   resid    = Y - X * reshape(aDraws, Nx, Ny);
end



function [X, Y] = lag4VAR(Y, p, constFlag)
% [X, Y] = lag4VAR(Y, p, constFlag)

% Elmar Mertens
% www.elmarmertens.ch

if nargin < 3
   constFlag = true;
end

[nobs, N] = size(Y);
T = nobs - p;
k = N * p;

if constFlag
   k = N * p + 1;
end

% construct regressors
X = ones(T, k);
for i = 1 : p
   X(:, (i - 1) * N + (1 : N))   = Y((p+1 : end) - i,       :);
end
if nargout > 1
   Y = Y(p+1 : end,:);
end
