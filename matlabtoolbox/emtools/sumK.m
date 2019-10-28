function Sk = sumK(x, horizons, exante)
% sumK rolling sum over matrix columns
% USAGE: function Sk = sumK(x, horizons, exante)
% this function allows for x being a matrix
%
% The output Sk has dimensions T x N x h where T x N are the dimensions of x
% and h is the number of elements in the vector horizons
%
% optional argument exante [=false]: stores the cumulated sum at the *beginning* of the summation period in Sk

%   Coded by  Elmar Mertens, em@elmarmertens.com

error(nargchk(2,3,nargin))
if nargin < 3
   exante = false;
end

horizons = horizons(:);

[T, N] = size(x);

Sk = NaN([T, N, length(horizons)]);

for i = 1 : length(horizons)
    k = horizons(i);
    if exante
       offset = k - 1;
    else
       offset = 0;
    end
    for t = 1 + (k - 1) : T 
        Sk(t-offset, :, i)  = sum(x(t - (k - 1) : t,:));
    end
end

Sk = squeeze(Sk);
