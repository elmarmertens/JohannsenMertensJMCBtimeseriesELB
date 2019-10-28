function [fstat, fprb, chi2stat, chi2prb] = linwald(resultu, R, r, EViewsflag, printflag)
% PURPOSE: computes Wald F-test for two regressions
%---------------------------------------------------
% USAGE: [fstat, fprb, chi2stat, chi2prb] = linwald(resultu, R, r, EViewsflag)
%    or: [fstat, fprb, chi2stat, chi2prb] = linwald(resultu) which calculates
%        standard F-test (assuming first regressor is constant)!
%  Note: if called without output arguments, it prints to the screen
% Where: resultu    = results structure from ols() unrestricted regression
%                     (accepts also hwhite etc., note: requires ols() with
%                     results structure yielding V matrix of estimates (as updated by EMT)
%        R          = full rank matrix of linear constraints R*b = r (q x nvar).
%                     R can be also a vector of indices into beta (or a logical index).
%                     Identifying which values of beta should be jointly tested against zero.
%        r          = vector of linear constraints coefficients R*b = r (q x nvar) [r = 0]
%        EViewsflag = Flag, if set to one calculates chi2 stat with dof-adjustment
%                     (i.e. chi2 = F/q) as EViews does [0].
%                     Note: To use EViewsflag while using default settings of R,r
%                     call linwald(resultu, R, [], EViewsflag) resp.
%                     linwald(resultu, [], [], EViewsflag)
%---------------------------------------------------
% RETURNS: fstat = {(T-K)/q * (R*b-r)' [R*(X'X)(^-1)*R']^-1 (R*b-r) / u'u}
%                  (see Hayashi or Amemiya 1985, p. 29, latter uses Q'b=c instead of R*b=r)
%          fprb  = marginal probability for fstat
% NOTE:  large fstat => reject the restrictions as inconsisent
%                       with the data
%---------------------------------------------------
% SEE ALSO: ols(), waldf(),
%---------------------------------------------------

%   Coded by  Elmar Mertens, em@elmarmertens.com

% written by:
% Elmar Mertens
% Study Center Gerzensee
% email: elmar.mertens@szgerzensee.ch

narginchk(1,5)
switch nargin
    case 1
        EViewsflag = 0;
        R = [];
        r = [];
        printflag = [];
    case 2
        EViewsflag = 0;
        r = [];
        printflag = [];
    case 3
        EViewsflag = 0;
        printflag = [];
    case 4
        printflag = [];
        %do nothing, all ok. I do not check validity
end

if isempty(printflag)
    printflag = 0;
end
if nargout == 0 % override any settings if otherwise no output
    printflag = 1;
end

if isempty(R)
    R = [zeros(resultu.nvar - 1, 1) eye(resultu.nvar - 1)];
else
    % Check whether R is vector and transform into matrix, if necessary
    if isvector(R)
        R = R(:);
        if islogical(R)
            R = find(R);
        end
        nr = rows(R);
        dummy = zeros(nr, resultu.nvar);
        for i = 1 : nr
            dummy(i, R(i)) = 1;
        end
        R = dummy;
    else
        nr = rows(R);
    end
end


if isempty(r)
    r = zeros(nr, 1);
end


if (cols(r) ~= 1) || (rows(r) ~= nr)
    error('linwald: Inconsistent size of arguments R, c')
end
if isstruct(resultu) == 0
    error('linwald requires an ols results structure as input');
end

%q equals number of restrictions
T   = resultu.nobs;
K   = resultu.nvar;
b   = resultu.beta;

if isfield(resultu, 'Tscale')
    if  resultu.Tscale        % all calculations below will assume that V was calculated without dof-adjustment
        V   = resultu.Vbeta;
    else
        V   = resultu.Vbeta * (T-K) / T;
    end
else
    V   = resultu.Vbeta;
end

% chi2stat = (R * b - r)' *  inv(R * V * R') * (R * b - r); %dimension are pretty small, hence Gauss-Jordian inverse won't harm
chi2stat = (R * b - r)'  / (R * V * R') * (R * b - r); 
fstat    = chi2stat  / nr * (T - K) / T;  % F-test will *here* always be calculated with dof-adjustment
fprb     = fdis_prb(fstat,nr, T - K);
if EViewsflag
    chi2stat = chi2stat * (T-K) / T;
end
chi2prb  = 1 - chis_prb(chi2stat, nr);


if printflag == 1
    fprintf(1,'\n');
    fprintf(1,'Wald Statistics\n');
    fprintf(1,'- F-statistic     = %16.8f \n',fstat);
    fprintf(1,'  probability     = %16.4f \n',fprb);
    fprintf(1,'  num,denom dof   = %4d,%4d\n',nr, T - K);
    fprintf(1,'- Chi2-statistic  = %16.8f \n',chi2stat);
    fprintf(1,'  probability     = %16.4f \n',chi2prb);
    fprintf(1,'  dof             = %4d \n', nr);
end
