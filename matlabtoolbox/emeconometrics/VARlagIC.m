function [AIC, SIC, HQIC, avgllf, Q, Qpvalue, Nsims, maxroots, TlessK, QpvalueChi2] = VARlagIC(Y, P, S, NQsims)
% function [AIC,SIC, HQIC, avgllf, Q, Qpvalue, Nsims, maxroots, TlessK, QpvalueChi2] = VARlagIC(Y, P, S, NQsims)
% reports information criteria for VAR over data Y with over grid "P" lags
% S are lags to be considered in LjungBox Portmanteau tests

%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 4 || isempty(NQsims)
   QsimFlag = false;
else
   % legacy check
   if islogical(NQsims)
      QsimFlag = NQsims;
      NQsims   = [];
   else
      QsimFlag = true;
   end
end

if nargout > 4
   if isscalar(S)
      S = 1 : S;
   end
end

if isscalar(P)
   P = 1 : P;
end


[AIC, SIC, HQIC, avgllf, maxroots, TlessK]   = deal(NaN(length(P), 1));
if nargout > 4
   [Q, Qpvalue, QpvalueChi2]                 = deal(NaN(length(P), length(S)));
end

for p = 1 : length(P)
   VAR            = VARlsdev(Y, P(p), [], true);
   TlessK(p)      = VAR.T - VAR.k;
   if TlessK(p) < 10
      warning('em:msg', 'you are looking at VARs with less than 10 dof!')
   end
   avgllf(p)      = VAR.llf / VAR.T;
   AIC(p)         = VAR.AIC;
   SIC(p)         = VAR.SIC;
   HQIC(p)        = VAR.HQIC;
   if nargout > 4
      [jim, jack]                = LjungBoxPortmanteau(VAR.e, S, P(p));
      Q(p, :)                    = jim';
      QpvalueChi2(p, :)          = jack';
      if QsimFlag
         percy          = [.9 .95 .99];
         Qsim           = VARmontecarloPortmanteau(VAR, NQsims, S);
         Nsims          = size(Qsim, 1);
         john           = repmat(Q(p,:), 3, 1) > quantile(Qsim, percy);
         Qpvalue(p,:)   = 1 - max(john .* repmat(percy', 1, length(S)));
      else
         Qpvalue(p,:)   = QpvalueChi2(p,:);
         Nsims = NaN;
      end
      
      if nargout > 7
         maxroots(p)     = max(abs(VARroots(VAR)));
      end
   end
   
end

