function [mu, tstat, pvalue, se] = dmtest(loss1, loss2, nlag)
% DMTEST ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 28-Jul-2017 20:32:48 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.2.0.556344 (R2017a) 
% FILENAME  : dmtest.m 


nanny  = isnan(loss1) | isnan(loss2);

delta = loss1(~nanny) - loss2(~nanny);
Nobs  = length(delta);

reggae = nwest(delta, ones(Nobs,1), nlag);

mu      = reggae.beta;
tstat   = reggae.tstat;
se      = sqrt(reggae.Vbeta);
pvalue  = tdis_prb(tstat,Nobs-1);

