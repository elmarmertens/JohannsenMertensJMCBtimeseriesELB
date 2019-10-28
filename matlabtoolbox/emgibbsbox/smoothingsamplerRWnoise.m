function [taudraws, tau0draw, noisedraws, z, z0, e] = smoothingsamplerRWnoise(y,vartrend,varnoise,tau00,V00,rndStream)
% SMOOTHINGSAMPLERRWNOISE ... 
%  
%   ... 

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 23-Nov-2011 11:35:01 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : smoothingsamplerRWNoise.m 

if nargin < 6
    rndStream = getDefaultStream;
end

y  = y(:);
T 	= length(y);

if isscalar(vartrend)
   tmp         = vartrend;
   vartrend    = zeros(T,1);
   vartrend(:) = tmp;
else
   vartrend = vartrend(:);
end
if isscalar(varnoise)
   tmp         = varnoise;
   varnoise    = zeros(T,1);
   varnoise(:) = tmp;
else
   varnoise = varnoise(:);
end
%% collect all randomdraws
z0 = randn(rndStream);
z  = randn(rndStream,T,1);
e  = randn(rndStream,T,1);

%% simulate plus
tau0plus = tau00 + sqrt(V00) * z0;
tauplus  = tau0plus + cumsum(z .* sqrt(vartrend));
yplus    = tauplus + e .* sqrt(varnoise);

ytilde   = yplus - y;


%% forward filter
[Sigma, SigmaStar]   = deal(NaN(T,1));
tautT                = NaN(T,1);

% initial distribution adjusted for ytilde
SigmaPrevious = V00; % still V00, in order to have identical gains
tauPrevious   = 0;   % this way, initial terms cancel out   

for t = 1 : T
   
   Sigma(t)       = SigmaPrevious + vartrend(t);
   gain           = Sigma(t) / (Sigma(t) + varnoise(t));
   SigmaStar(t)   = (1 - gain)^2 * Sigma(t) + gain^2 * varnoise(t);
   
   tautT(t)       = (1 - gain) * tauPrevious + gain * ytilde(t);
   
   % prepare next loop
   tauPrevious       = tautT(t); 
   SigmaPrevious     = SigmaStar(t); 
   
end


%% backward filter
for t = T - 1 : -1 : 1
   gain     = SigmaStar(t) / Sigma(t+1);
   tautT(t) = (1 - gain) * tautT(t) + gain * tautT(t+1);
end
gain  = V00 / Sigma(1);
tau0T = gain * tautT(1); % again: treat tau00 as zero

%% put draws together
taudraws    = tauplus - tautT;
if nargout > 1
    tau0draw    = tau0plus - tau0T;
    if nargout > 2
        noisedraws  = y - taudraws;
    end
end
