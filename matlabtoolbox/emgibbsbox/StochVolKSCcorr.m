function [h, h0, kai2States] = StochVolKSCcorr(logy2, h, hVCV, Eh0, Vh0, KSC, KSCt, Nsv, T, rndStream)
% StochVolKSC performs a Gibbs updating step on a SV model and it works over 
% Nsv independent SV residuals
%
% Uses Kim, Shephard and Chib normal mixtures
%
% USAGE: [h, kai2States] =  StochVolKSCabc(logy2, kai2States, hInno, Eh0, Vh0, 
%                           KSC, KSCt, Nsv, T, rndStream)
%
% See also abcDisturbanceSmoothingSampler, getKSCvalues

%   Coded by  Elmar Mertens, em@elmarmertens.com


%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 28-Aug-2009 12:07:01 $
% $Revision : 1.00 $
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : StochVolKSC.m.m

if isscalar(Eh0)
    Eh0 = repmat(Eh0, Nsv, 1);
end
if isscalar(Vh0)
    Vh0 = repmat(Vh0, Nsv, 1);
end

%% CORRIGENDUM CHANGES ORDER OF GIBBS STEPS!

%% draw mixture states
% zdraws are standardized draws for each component of the normal mixture 
% zdraws is thus Nsv x T x 7 
zdraws      = bsxfun(@minus, logy2 - h, KSCt.mean) ./ KSCt.vol;

% construct CDF
% factor of sqrt(2 * pi) can be ommitted for kernel
pdfKernel           = KSCt.pdf ./ KSCt.vol .* exp(-.5 * zdraws.^2); 
cdf                 = cumsum(pdfKernel, 3);                % integrate
cdf(:,:,1:end-1)    = bsxfun(@rdivide, cdf(:,:,1:end-1), cdf(:,:, end)); 
% switch size(cdf, 3) 
%    case 7
%       cdf(:,:,1:end-1)    = cdf(:,:,1:end-1) ./ cdf(:,:, [7 7 7 7 7 7]);   % normalize
%    case 2
%       cdf(:,:,1)    = cdf(:,:,1) ./ cdf(:,:,2);   % normalize
%    otherwise
%       error('KSC error: currently, only 7 or 2 components implemented')
% end
cdf(:,:,end)        = 1;    % normalize

% bsxfun appears to be pretty slow
% check       = bsxfun(@rdivide, cdf, cdf(:,:,end));  % normalize
% checkdiff(check, cdf);

% draw states
kai2States  = sum(bsxfun(@gt, rand(rndStream, Nsv, T), cdf), 3) + 1;
obs         = logy2 - KSC.mean(kai2States);

%% KSC State Space
A     = eye(Nsv);
B     = chol(hVCV)';
C     = eye(Nsv);
sqrtR = zeros(Nsv,Nsv,T);
for n = 1 : Nsv
    sqrtR(n,n,:) = KSC.vol(kai2States(n,:));
end

sqrtVh0 = diag(sqrt(Vh0));

% [h, hshock, h0] = abcDisturbanceSmoothingSampler(A, B, C, obs, Eh0, sqrtVh0, 1, ...
%     sqrtR, [], rndStream); %#ok<ASGLU>
[h, hshock, h0] = abcrDisturbanceSmoothingSampler1draw(A, B, C, obs, Eh0, sqrtVh0, ...
    sqrtR, rndStream); %#ok<ASGLU>



