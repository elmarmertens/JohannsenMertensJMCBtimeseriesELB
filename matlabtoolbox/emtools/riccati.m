function [X, F, converged, message] = riccati(A,B,Q,R,N,X0,tolerance,maxiter)
% RICCATI solves riccati equation via value function iteration
% similar to matlab's dare, riccati solves
%
% X = Q + A' X A - (A' X B + N) * (B' X B + R)^(-1) * (B' X A + N')
%
% the gain matrix is F = -(R + B' * X * B) \ (N' + B' * X * A)
%
% the only difference with matlab's dare is that matlab notation says S 
% instead of N and that matlab's dare allows for a sixth input matrix E
%
% [X, F, converged, message] = riccati(A,B,Q,R,N,X0,tolerance,maxiter)
%
% where tolerance and maxiter are optional with defaults of 1e-12 and 1e4
% and converged is a boolean variable equal to true if the VFI has
% successfully converged
%
% When RICCATI is called with less than two output arguments, an error
% message is generated in case of non-convergence
%
% The algorithm is identical to the one used by linregulatorVFI
%
% see also dare, linregulatorVFI

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 24-Nov-2010 15:07:03 $
% $Revision : 1.00 $
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : riccati.m


if nargin < 4 || isempty(R)
    R = zeros(size(B, 2));
end
if nargin < 5 || isempty(N)
    N = zeros(size(A, 2), size(B, 2));
end

if nargin < 6 || isempty(X0)
    X0 = Q; % zeros(size(A));
end

if nargin < 7 || isempty(tolerance)
    tolerance	= 1e-12;
end

if nargin < 8 || isempty(maxiter)
    maxiter     = 1e6;
end


% if nargout > 4
%    history = repmat(struct('F', NaN(size(B, 2), size(B, 1)), 'X', repmat(NaN, size(A))), maxiter, 1);
% else
%    history = [];
% end

iter        = 0;
converged   = false;

while ~converged && iter < maxiter
    
    iter = iter + 1;
    
    F  = -(R + B' * X0 * B) \ (N' + B' * X0 * A);
    AA = A + B * F;
    X  = Q + N * F + F' * N' + F' * R * F + AA' * X0 * AA;
    
    dV = X - X0;
    converged = max(abs(dV(:))) < tolerance;
    
    X0 = X;
    
    %    if ~isempty(history)
    %       history(iter).F = F;
    %       history(iter).X = X;
    %    end
    
end

% if ~isempty(history)
%    history = history(1:iter);
% end

if ~converged
    message = sprintf('RICCATI failed to converge after %d iteration, mae is %f', iter, max(abs(dV(:))));
    if nargout < 3
        error(message) %#ok<SPERR>
    end
else
    message = [];
end

