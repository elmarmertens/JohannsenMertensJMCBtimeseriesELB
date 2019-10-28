function haccme = haccme(u,nlag)
% PURPOSE: computes VCV adjusted for autocorrelation following Newey and West(1987)
%          The case of general heteroskedasticity without serial correlation 
%          of White (1980) is obtained for nlag=0 .
%----------------------------------------------------------------------------------
% USAGE: haccme       = HACCME(u,nlag)
% where: u            = north x nobs sample moment residuals
%                       also knwon as hhat (Newey and West) or f (Hansen and Singleton)
%                       i.e. u(:,t) is an observation of sample moment function
%        nlag         = scalar, denotes the lags to use (0 for White HCCME)
%----------------------------------------------------------------------------------
% RETURNS: haccme = north x north HAC VCV (HACCME) of the sample moments g=mean(u,2)
%----------------------------------------------------------------------------------
% References: Newey and West, 1987, Econometrica
%             also H. White 1980, Econometrica Vol. 48 pp. 818-838. 
%             Davidson and MacKinnon, 1993, pp. 553f. 
%----------------------------------------------------------------------------------

%   Coded by  Elmar Mertens, em@elmarmertens.com

% written by
% Elmar Mertens
% Graduate Student at the University of St. Gallen
% Hoehenweg 13
% CH-9000 St. Gallen
% Switzerland
% +41 (0)71 223 8961
% elmar.mertens@student.unisg.ch% 

% based on mcov.m written by
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu

% ---------------------------------------------------------------------------------
% CHECKING INPUT AND INTITALIZING VARIABLES
% ---------------------------------------------------------------------------------

error(nargchk(1,2,nargin));
if nargin<2
   nlag=0;
end

if nlag<0
   warning('em:warn','nlag is negative - using abs(nlag) instead');
   nlag=abs(nlag);
end
% [north nobs]=size(u);
nobs = size(u, 2);

% if nlag > 10
%    warning('em:msg', 'don''t you want to make a dof adjustment?')
% end

% ---------------------------------------------------------------------------------
% THE CALCULATION STARTS HERE
% ---------------------------------------------------------------------------------
haccme= u * u'; % that's omega_zero
for j = 1 : nlag;
	% omega_j=zeros(north,north);
   w  = 1 - j / (nlag + 1);
%    u1 = u(:,(j + 1) : nobs);
%    u2 = u(:, 1 : nobs - j);
%    omega_j = u1 * u2';
   omega_j = u(:,(j + 1) : nobs) *  u(:, 1 : nobs - j)';
   haccme = haccme + w * (omega_j + omega_j');
end
haccme = haccme / nobs; 

% Concerning the candidate finite sample correction n / (N-K) DMK, p.610, remark that this topic has not much been investigated "at the timw of writing"
