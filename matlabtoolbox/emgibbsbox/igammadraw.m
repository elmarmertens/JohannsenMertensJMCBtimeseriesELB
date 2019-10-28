function [igams, gams] = igammadraw(sigmaT, dof, Ndraw, rndStream)
% IGAMMADRAW produces random number from (multiple) inverse gamma distributions
%
% USAGE: [igams, gams] = igammadraw(sigmaT, dof, Ndraw, rndStream)
% NOTATION: mean of inverse Gamma is sigmaT / (dof - 2)
% output is length(sigmaT) x Ndraw
% 
% NOTE: THIS IS RATHER A UNIVARIATE INVERSE WISHART (and does not follow the notation for the inverse gamma)
% the proper maping from the univaraite-inverse Wishart to the igamma would
% be alpha = dof /2 and beta = sigmaT / 2
% See also: iwishdraw

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 04-Mar-2009 15:56:44 $
% $Revision : 1.00 $
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : igammadraw.m

%% parse inputs
if nargin < 2
    dof = [];
end
if nargin < 3
    Ndraw = 1;
end
if nargin < 4
    rndStream = getDefaultStream;
end

if ~isscalar(dof)
    if length(sigmaT) ~= length(dof)
        error('dof is not a scalar, but has different size than sigmaT')
    end
    igams = NaN(length(dof),Ndraw);
    if nargout == 1
        for n = 1 : length(dof)
            igams(n,:) = igammadraw(sigmaT(n), dof(n), Ndraw, rndStream);
        end
    else
        gams = NaN(length(dof),Ndraw);
        for n = 1 : length(dof)
            [igams(n,:), gams(n,:)] = igammadraw(sigmaT(n), dof(n), Ndraw, rndStream);
        end
    end
else
    
    %% check dof
    sigmaT  = sigmaT(:);
    N       = length(sigmaT);
    
    if isempty(dof)
        dof = 1;
    end
    
    if dof < 1
        warning('em:msg', 'dof must be larger than 0; will be set to 1. (dof = %d)', dof)
        dof = 1;
    end
    
    %% draw random normal numbers
    z  = randn(rndStream, N, Ndraw, dof);
    
    %% construct inverse-gammas
    zz       = squeeze(sum(z.^2,3));
    igams    = bsxfun(@rdivide, sigmaT, zz);
    
    %% invert to gammas (if output is required)
    if nargout > 1
        gams = 1 ./ igams;
    end
end

