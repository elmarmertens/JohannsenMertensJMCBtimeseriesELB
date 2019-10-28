function X = dlyapdoubling(A,B,Q,tolerance,maxiter)
% DLYAPDOUBLING solves the discrete-time Lyapunov function with a doubling algorithm
% this code mimicks Matlab's own dlyap *for pure Lyapunov equations*
% 
% X = dlyapdoubling(A,Q) solves X = A * X * A' + Q
%
% X = dlyapdoubling(A,B,Q) solves the Sylvester equation X = A * X * B' + Q 
%
% Optional parameters: tolerance and maxiter
% which can be specified via dlyap(A,Q,tolerance) and dlyap(A,Q,tol,maxiter) 
% (or dlyap(A,B,Q,tolerance) and dlyap(A,B,Q,tolerance,maxiter)). 
% The default values are tolerance = 1e-12 and maxiter = 1e4
%
% See also dlyapvec

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 24-Nov-2010 11:58:24 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : dlyap.m 

switch nargin
    case 2 % dlyap(A,B)
        Q = B;
        B = A;
        tolerance = [];
        maxiter   = [];
    case 3 % dlyap(A,B,tolerance) or dlyap(A,B,Q)
        if isscalar(Q) 
            tolerance   = Q;
            Q           = B;
        else
            tolerance = [];
        end
        maxiter   = [];
    case 4 % dlyap(A,B,tolerance,maxiter) or % dlyap(A,B,Q,tolerance)
        if isscalar(Q) 
            maxiter     = tolerance;
            tolerance   = Q;
            Q           = B;
        end
end

if isempty(tolerance)
   tolerance = 1e-12;
end
if isempty(maxiter)
   maxiter = 1e4;
end

X0 = Q;
AA = A;
BB = B;

converged   = false;
iter        = 0;

while ~converged && iter < maxiter 
    
    iter = iter + 1;
    
    X  = X0 + AA * X0 * BB';
    AA = AA * AA;           % this is the doubling step
    BB = BB * BB;           % this is the doubling step

    dX          = X - X0;
    converged   = max(abs(dX(:))) < tolerance; 
    X0          = X;

    if any(isnan(X0(:)))
        break
    end
    
end

if ~converged
    error('dlyapdoubling failed to converge after %d iterations.\n\t largest dX was %e\n\t tolerance %e\n\nLargest eigenvalue of transition is %e\n\n', iter, max(abs(dX(:))), tolerance, max(abs(eig(A))))
end

% X = (X + X') / 2; % enforce symmetry
