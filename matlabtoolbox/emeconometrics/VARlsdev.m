function results = VARlsdev(Y, p, dfFlag, quickFlag)
% Efficient Computation of VARE -- demeaning variables first
%---------------------------------------------------
% USAGE:    results = VARlsdev(y, p, dfFlag, quickFlag))
% where:    y       = an (T x N) matrix of y-vectors
%           p       = the lag length
%---------------------------------------------------
% 
% This function always includes a constant in the VAR. 
% The inverse of X' * X is computed from partioned matrix formulas, 
% which avoids ill-conditioning due to large mean values of the data
%
% See also VARls

%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 3 || isempty(dfFlag)
   dfFlag = false;
end
if nargin < 4 || isempty(quickFlag)
   quickFlag = false; % EViews 5 is now compatible with this, for EViews 4, set to true
end

results.meth      = 'VARls';
results.constFlag = true;    % always do constant, this entry accomodates legacy code

[nobs, N] = size(Y);
T = nobs - p;
k = N * p + 1;  

% construct regressors
if ~quickFlag
   [X, results.CompanionState] = deal(ones(T, k));
   for i = 1 : p;
      results.CompanionState(:, (i - 1) * N + (1 : N))   = Y((p+1 : end) - i + 1,   :);
      X(:, (i - 1) * N + (1 : N))                        = Y((p+1 : end) - i,       :);
   end
else
   X = ones(T, k);
   for i = 1 : p;
      X(:, (i - 1) * N + (1 : N))               = Y((p+1 : end) - i,       :);
   end
end

Y0    = Y(1 : p, :);
Y     = Y(1 + p : end, :);
Xbar  = mean(X(:,1:end-1))';

% Xdev  = X(:,1:end-1) - ones(T,1) * Xbar';
Xdev  = bsxfun(@minus, X(:,1:end-1), Xbar');

XXdev       = Xdev' * Xdev;


if quickFlag
   XXidev   = XXdev \ eye(k - 1); % for small data (ca 100) faster by factor 2-3.5 than qr, for ca 1000 factor 1.5-2, see JPL manual. Accuracy is worse  for ill-conditioned problems though
else
   condi = rcond(XXdev);
   if condi < sqrt(eps)
      warning('em:VARls', 'Badly Conditioned System: rcond(XX) = %e', condi)
   end
   [tmp r]  = qr(Xdev,0);
   XXidev   = (r'*r) \ eye(k - 1);
end


XXidevXbar  = -XXidev * Xbar;
XXi         = [XXidev XXidevXbar; XXidevXbar' (1/T + Xbar' * XXidev * Xbar)];

% XXi2  = [XXidev -XXidev * Xbar; - Xbar' * XXidev (1/T + Xbar' * XXidev * Xbar)];
% checkdiff(XXi, XXi2);
% checkdiff(XXi, inv(X' * X));

PI    = XXi * X' * Y; % compute coefficients
e     = Y - X * PI;

if dfFlag
   Omega = (e' * e) / (T - k);
else
   Omega = (e' * e) / T;
end
if ~quickFlag
   aVar  = kron(Omega, XXi * T);
   if dfFlag
      SE    = sqrt(diag(aVar) / (T - k));
   else
      SE    = sqrt(diag(aVar) / T);
   end
   SE = reshape(SE, size(PI));
else
   aVar = NaN;
   SE   = NaN;
end

results.k = k;
results.K = k * N;
results.T = T;
results.N = N;
results.p = p;

results.e               = e; 
results.PI              = PI;
results.X               = X;
results.Y               = Y;
results.Y0              = Y0;

% redundant:
% results.data = Y; 

results.Omega     = Omega;
results.aVar      = aVar;
results.SE        = SE;

% companion
results.F       = [PI(1 : N * p, :)'; eye(N*(p-1), N*p)];
results.F0      = [PI(end, :)'; zeros(N*(p-1), 1)];
if ~quickFlag
   results.roots   = eig(results.F);
end
results.G       = eye(N * p, N);
% results.Omega   = G * Omega * G';

results.llf = - (T * N / 2) * (1 + log(2 * pi)) - T / 2 * log(det(Omega));
results.dfFlag = dfFlag;

[results.AIC, results.SIC, results.HQIC] = IC(results.llf, results.T, results.K);
