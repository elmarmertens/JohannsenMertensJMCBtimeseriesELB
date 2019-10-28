function [iwishs, cholSigmaT, wishs] = iwishdraw(SigmaT, dof, Ndraw, cholSigmaT, rndStream)
% IWISHDRAW produces random number from inverse wishart distribution
% 
% USAGE: [iwishs, cholSigmaT, wishs] = iwishdraw(SigmaT, dof, Ndraw, cholSigmaT, rndStream)
%
% NOTATION: mean of inverse Wishart is SigmaT / (dof - N - 1)), where N = length(Sigma))
%
% Notice, for a proper draw (with first moment defined) , need dof > N + 1. 
%   ...

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 04-Mar-2009 15:56:44 $
% $Revision : 1.00 $
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : iwishdraw.m

%% parse inputs
if nargin < 2 
    dof = [];
end
if nargin < 3
   Ndraw = 1;
end
if nargin < 4
   cholSigmaT = [];
end
if nargin < 5
   rndStream = getDefaultStream;
end

%% compute choleski (if necessary)
if isempty(cholSigmaT)
    cholSigmaT = chol(SigmaT)';
end

%% check dof
N = size(cholSigmaT, 1);
if isempty(dof)
   dof = N;
end
% if dof < N
%    %    warning('em:msg', 'dof must be larger or equal than N, will be set to N. (dof = %d, N = %d)', dof, N)
%    dof = N;
% end

%% allocate memory for results
iwishs  = zeros(N, N, Ndraw);
if nargout > 2
   wishs       = zeros(N, N, Ndraw);
   icholSigmaT = eye(N) / cholSigmaT;
end

%% draw random normal numbers
z       = randn(rndStream, N, dof, Ndraw);

%% loop over Ndraw and construct inverse-wisharts
if Ndraw == 1
    zz      = z * z';
    iwishs  = cholSigmaT / zz * cholSigmaT';
    if nargout > 2
        wishs	= icholSigmaT' * zz * icholSigmaT;
    end
else
    for n = 1 : Ndraw
        zz              = z(:,:,n) * z(:,:,n)';
        iwishs(:,:,n)   = cholSigmaT / zz * cholSigmaT';
        if nargout > 2
            wishs(:,:,n) = icholSigmaT' * zz * icholSigmaT; %  checkdiff(wishs(:,:,n), inv(iwishs(:,:,n)));
        end
    end
end
