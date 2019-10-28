function [flag, delta] = checkdiff(a,b,tol,msg)
% function [flag, delta] = checkdiff(a,b,tol,msg)
% flag = true if there is a significant difference (in mean absolute error sense) between a and b
% tol is significance limit (default tol = 1e-10)
%
% if nargout = 0, checkdiff returns a warning message when there is a difference
% the message can be customized using the msg argument
% when msg = [], checkdiff tries to recognize a and b via inputname(1) and inputname(2)

%   Coded by  Elmar Mertens, em@elmarmertens.com

% Elmar Mertens
% www.elmarmertens.ch

% error(nargchk(2,4,nargin))
alabel = inputname(1);
if nargin < 2 || isempty(b);
   b = zeros(size(a));
   blabel = '';
else
   blabel = inputname(2);
end


if nargin < 3 || isempty(tol)
   tol = 1e-10;
end
if nargin < 4
   msg = [];
end


if ~isnumeric(tol)
   error('tol must be numeric') %#ok<ERTAG>
end

if isnumeric(a) && isnumeric(b)
   
   a           = a(:);
   b           = b(:);
   nanny       = ~any(isnan([a b]), 2);
   delta       = max(abs(a(nanny) - b(nanny)));
   flag        = delta > tol;
   
else
   
   flag      = ~isequal(a,b);
   delta     = NaN;
   
end

if isempty(flag) % there have been only NaN's
   flag = false;
end

if flag && (nargout == 0 || ~isempty(msg))
   if isempty(msg)
      warning('em:checkdiff', 'mae(%s,%s)=%e', alabel, blabel, delta)
   else
      warning('em:checkdiff', '%s -- mae=%e', msg, delta)
   end
end
