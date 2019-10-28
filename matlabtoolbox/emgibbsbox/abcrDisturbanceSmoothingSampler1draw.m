function [Xdraws, disturbanceDraws, X0draws, noiseDraws] = ...
    abcrDisturbanceSmoothingSampler1draw(A, B, C, Ydata, X00, cholSigma00, sqrtR, rndStream)
% ABCDISTURBANCESMOOTHINGSAMPLER
% ....

% accepts only 2D: ABC; only sqrtR can be 3D 

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 08-Aug-2009 17:58:16 $
% $Revision : 1.00 $
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : abcDisturbanceSmoothingSampler.m


%% parse inputs
Nx                = size(A, 1);
[Ny, T]           = size(Ydata);
Nw                = size(B,2);

% if nargin < 8
%     sqrtR = [];
% end

% if nargin < 9
%     rndStream = getDefaultStream;
% end

%% init Variables and allocate memory
% if ~isempty(sqrtR) && ismatrix(sqrtR)
%     sqrtR = repmat(sqrtR, [1 1 T]);
% end

I                 = eye(Nx);
Iy                = eye(Ny);

%% allocate memory
[Sigmattm1, ImKC]           = deal(zeros(Nx, Nx, T));
invSigmaYttm1               = zeros(Ny, Ny, T);
Ytilde                      = zeros(Ny, T);
[XtT, Xttm1, Xplus]         = deal(zeros(Nx, T));


%% generate plus data

wplus   = randn(rndStream, Nw, T);
eplus            = randn(rndStream, Ny, T);
X0plus = X00 + cholSigma00 * randn(rndStream, Nx, 1); 

%% Forward Loop: Kalman Forecasts
[Sigma00, Sigmatt] = deal(cholSigma00 * cholSigma00');
Xtt     = zeros(Nx,1); % use zeros, since projection on difference between Y and Yplus
% BSigmaB = zeros(Nx, Nx, T);

disturbanceplus  = zeros(Nx, T);
noiseplus        = zeros(Ny, T);

BSigmaB   = B * B';

for t = 1 : T
    
    % "plus" States and priors
    disturbanceplus(:,t)  = B * wplus(:,t);
    
    if t == 1
        Xplus(:,t) = A * X0plus + disturbanceplus(:,t);
    else
        Xplus(:,t) = A * Xplus(:,t-1) + disturbanceplus(:,t);
    end
    
    % priors
    Sigmattm1(:,:,t)        = A * Sigmatt * A' + BSigmaB;
    Xttm1(:,t)              = A * Xtt;
    
    
    
    % observed innovation
    noiseplus(:,t)    = sqrtR(:,:,t) * eplus(:,t); 
    Yplus             = C * Xplus(:,t) + noiseplus(:,t);
    SigmaYttm1        = C * Sigmattm1(:,:,t) * C' + ...
        sqrtR(:,:,t) * sqrtR(:,:,t)';
        
   
    Ytilde(:,t) = Ydata(:,t)  - Yplus  - C * Xttm1(:,t);
    
    % Block Inverse of Y-VCV, accounting for missing obs (with zero VCV)
    invSigmaYttm1(:,:,t) = Iy / SigmaYttm1;

    % Kalman Gain
    K                       = (Sigmattm1(:,:,t) * C') * invSigmaYttm1(:,:,t);
    ImKC(:,:,t)             = I - K * C;
    
    % posteriors
    Sigmatt                 = ImKC(:,:,t) * Sigmattm1(:,:,t);
    
    Xtt                     = Xttm1(:,t) + K * Ytilde(:,t);
   
end

%% Backward Loop: Disturbance Smoother
XtT(:,T)        = Xtt;

StT             = C' * (invSigmaYttm1(:,:,T) * Ytilde(:,T));

if nargout > 1
    disturbancetT               = zeros(Nx, T);
    disturbancetT(:,T)          = BSigmaB * StT;
else
    disturbancetT        = [];
end

if nargout > 2 
    noisetT            = zeros(Ny, T);
    noisetT(:,T)       = Ytilde(:,T) - C * (XtT(:,T) - Xttm1(:,T));
else
    noisetT      = [];
end


for t = (T-1) : -1 : 1
    Atilde      = A * ImKC(:,:,t);
    
    StT         = Atilde' * StT + ...
        C' * (invSigmaYttm1(:,:,t) * Ytilde(:,t));
    XtT(:,t)    = Xttm1(:,t) + Sigmattm1(:,:,t) * StT;
    
    if ~isempty(disturbancetT)
        disturbancetT(:,t)        = BSigmaB * StT;
    end
    if ~isempty(noisetT)
        noisetT(:,t)       = Ytilde(:,t) - C * (XtT(:,t) - Xttm1(:,t));
    end
    
end

%% sample everything together (and reorder output dimensions)
Xdraws  = Xplus + XtT;

if nargout > 1
    
    disturbanceDraws  = disturbanceplus + disturbancetT;
    
    if nargout > 2
        
        X0T      = Sigma00 * A' * StT; % note: no mean added to X0T since it is already included in X0plus
        X0draws  = X0plus + X0T;
        
        
        if nargout > 3
            noiseDraws  = noiseplus + noisetT;
        end
    end
end


       
