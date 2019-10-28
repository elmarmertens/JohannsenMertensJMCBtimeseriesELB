function [AIC, SIC,HQIC] = IC(llf, T, K)
% function [AIC, SIC, HQIC] = IC(llf, T, K)

%   Coded by  Elmar Mertens, em@elmarmertens.com

% Elmar Mertens
% www.elmarmertens.ch

if (nargin == 1) && isstruct(llf)
    T   = llf.T;
    K   = llf.K;
    llf = llf.llf;
end

AIC   = - 2 * llf / T + 2 * K / T;
SIC   = - 2 * llf / T + K * log(T) / T;
HQIC  = - 2 * llf / T + 2 * K * log(log(T)) / T;
