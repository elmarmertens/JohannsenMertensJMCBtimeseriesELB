function [Xdraws, disturbanceDraws, X0draws, noiseDraws] = ...
    abcDisturbanceSmoothingSampler1draw(A, B, C, Ydata, X00, cholSigma00, sqrtR, rndStream)
% ABCDISTURBANCESMOOTHINGSAMPLER
% ....

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

if nargin < 7
    sqrtR = [];
end

if nargin < 8
    rndStream = getDefaultStream;
end

%% init Variables and allocate memory
if ~isempty(sqrtR) && ismatrix(sqrtR)
    sqrtR = repmat(sqrtR, [1 1 T]);
end


if ismatrix(A)
    A = repmat(A, [1 1 T]);
end
if ismatrix(B)
    B = repmat(B, [1 1 T]);
end
if ismatrix(C)
    C = repmat(C, [1 1 T]);
    % warning('em:msg', 'C should be three dimensional when there are missing data')
end

I                 = eye(Nx);
Iy                = eye(Ny);

%% allocate memory
[Sigmattm1, ImKC]           = deal(zeros(Nx, Nx, T));
invSigmaYttm1               = zeros(Ny, Ny, T);
Ytilde                      = zeros(Ny, T);
[XtT, Xttm1, Xplus]         = deal(zeros(Nx, T));


%% generate plus data

wplus   = randn(rndStream, Nw, T);
if ~isempty(sqrtR) 
    Nmeasurmentnoise = size(sqrtR, 2);
    eplus            = randn(rndStream, Nmeasurmentnoise, T);
end
X0plus = X00 + cholSigma00 * randn(rndStream, Nx, 1); 

%% Forward Loop: Kalman Forecasts
[Sigma00, Sigmatt] = deal(cholSigma00 * cholSigma00');
Xtt     = zeros(Nx,1); % use zeros, since projection on difference between Y and Yplus
BSigmaB = zeros(Nx, Nx, T);

disturbanceplus  = zeros(Nx, T);
if isempty(sqrtR)
    noiseplus        = [];
else
    noiseplus        = zeros(Ny, T);
end

for t = 1 : T
    
    % "plus" States and priors
    disturbanceplus(:,t)  = B(:,:,t) * wplus(:,t);
    BSigmaB(:,:,t)        = B(:,:,t) * B(:,:,t)';

    
    if t == 1
        Xplus(:,t) = A(:,:,t) * X0plus + disturbanceplus(:,t);
    else
        Xplus(:,t) = A(:,:,t) * Xplus(:,t-1) + disturbanceplus(:,t);
    end
    
    % priors
    Sigmattm1(:,:,t)        = A(:,:,t) * Sigmatt * A(:,:,t)' + BSigmaB(:,:,t);
    Xttm1(:,t)              = A(:,:,t) * Xtt;
    
    
    
    % observed innovation
    if isempty(sqrtR)
        
        Yplus            = C(:,:,t) * Xplus(:,t);
        SigmaYttm1       = C(:,:,t) * Sigmattm1(:,:,t) * C(:,:,t)';
        
    else
        
        noiseplus(:,t)    = sqrtR(:,:,t) * eplus(:,t); %#ok<AGROW>
        Yplus             = C(:,:,t) * Xplus(:,t) + noiseplus(:,t);
        SigmaYttm1        = C(:,:,t) * Sigmattm1(:,:,t) * C(:,:,t)' + ...
            sqrtR(:,:,t) * sqrtR(:,:,t)';
        
    end
   
    Ytilde(:,t) = Ydata(:,t)  - Yplus  - C(:,:,t) * Xttm1(:,t);
    
    % Block Inverse of Y-VCV, accounting for missing obs (with zero VCV)
    invSigmaYttm1(:,:,t) = Iy / SigmaYttm1;

    % Kalman Gain
    K                       = (Sigmattm1(:,:,t) * C(:,:,t)') * invSigmaYttm1(:,:,t);
    ImKC(:,:,t)             = I - K * C(:,:,t);
    
    % posteriors
    Sigmatt                 = ImKC(:,:,t) * Sigmattm1(:,:,t);
    
    Xtt                     = Xttm1(:,t) + K * Ytilde(:,t);
   
end

%% Backward Loop: Disturbance Smoother
XtT(:,T)        = Xtt;

StT             = C(:,:,T)' * (invSigmaYttm1(:,:,T) * Ytilde(:,T));

if nargout > 1
    disturbancetT               = zeros(Nx, T);
    disturbancetT(:,T)          = BSigmaB(:,:,T) * StT;
else
    disturbancetT        = [];
end

if nargout > 2 && ~isempty(sqrtR)
    noisetT            = zeros(Ny, T);
    noisetT(:,T)       = Ytilde(:,T) - C(:,:,T) * (XtT(:,T) - Xttm1(:,T));
else
    noisetT      = [];
end


for t = (T-1) : -1 : 1
    Atilde      = A(:,:,t+1) * ImKC(:,:,t);
    
    StT         = Atilde' * StT + ...
        C(:,:,t)' * (invSigmaYttm1(:,:,t) * Ytilde(:,t));
    XtT(:,t)    = Xttm1(:,t) + Sigmattm1(:,:,t) * StT;
    
    if ~isempty(disturbancetT)
        disturbancetT(:,t)        = BSigmaB(:,:,t) * StT;
    end
    if ~isempty(noisetT)
        noisetT(:,t)       = Ytilde(:,t) - C(:,:,t) * (XtT(:,t) - Xttm1(:,t));
    end
    
end

%% sample everything together (and reorder output dimensions)
Xdraws  = Xplus + XtT;

if nargout > 1
    
    disturbanceDraws  = disturbanceplus + disturbancetT;
    
    if nargout > 2
        
        X0T      = Sigma00 * A(:,:,1)' * StT; % note: no mean added to X0T since it is already included in X0plus
        X0draws  = X0plus + X0T;
        
        
        if nargout > 3
            noiseDraws  = noiseplus + noisetT;
        end
    end
end


       
