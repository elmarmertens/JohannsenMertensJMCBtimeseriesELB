function COR = cov2cor(COV, diagVol)
% function COR = cov2cor(COV)
% converts a covariance matrix into a correlation matrix
  
%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 2
    diagVol = false;
end

vols        = sqrt(diag(COV));
stochastic  = ~(vols == 0);
COR         = zeros(size(COV));

COR(stochastic,stochastic)  ...
            = COV(stochastic,stochastic) ./ (vols(stochastic) * vols(stochastic)');
if diagVol
    diagndx      = logical(eye(size(COV)));
    COR(diagndx) = vols;
end
        
        
