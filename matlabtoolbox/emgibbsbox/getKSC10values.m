function [KSC, KSCt] = getKSC10values(T, Nsv)
% GETKSCVALUES returns coefficients from Kim, Shephard and Chibs normal-mixture approximation of a chi2 variable
% KSC = getKSCvalues returns a structure with elements mean, vol, and pdf
% [KSC, KSCt] = getKSCvalues(T, Nsv) returns in addition, the structure KSCt, with the same fieldnames as KSC
% The fields of KSCt are "blown up" to dimension T x Nsv x 7 

%   Coded by  Elmar Mertens, em@elmarmertens.com


%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 28-Aug-2009 11:29:39 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : getKSCvalues.m 

if nargin < 2
    Nsv = 1;
end

KSC.mean    = [1.92677 1.34744 0.73504 0.02266 -0.85173 -1.97278 -3.46788 -5.55246 -8.68384 -14.65000];
KSC.var     = [0.11265 0.17788 0.26768 0.40611 0.62699 0.98583 1.57469 2.54498 4.16591 7.33342];
KSC.vol     = sqrt(KSC.var);
KSC.pdf     = [0.00609 0.04775 0.13057 0.20674 0.22715 0.18842 0.12047 0.05591 0.01575 0.00115];
KSC.cdf     = cumsum(KSC.pdf);

% blowup to cover time dimension
if nargout > 1
    fn = fieldnames(KSC);
    for f = 1 : length(fn)
        KSCt.(fn{f}) = repmat(permute(KSC.(fn{f}), [1 3 2]), [Nsv, T, 1]);
    end
end
