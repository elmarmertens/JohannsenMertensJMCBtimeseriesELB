%% preamble script for common settings to create plots for paper and appendix

%% load input data
datadir = pwd; % default: current directory; change to location where results are stored 

if doYGAP
    dataset    = 'spectreTB3MSGS020510OutputGapHeadline2018Q4';
    modellabel = 'baselineYgap';
else
    dataset    = 'spectreTB3MSGS020510Headline2018Q4';
    modellabel = 'baselineUgap';
end
Nyield = 3;

y       = importdata(fullfile(datadir, sprintf('%s.yData.txt', dataset)));
dates   = importdata(fullfile(datadir, sprintf('%s.dates.txt', dataset)));
yNaN    = logical(importdata(fullfile(datadir, sprintf('%s.yNaN.txt', dataset))));
[T, Ny] = size(y);

Ngap = Ny;

yLabel   = {'cgap','inflation', 'policy-rate', 'y2', 'y5', 'y10'};
ynames   = {'\tilde\pi_t', '\tilde c_t', '\tilde s_t', '\tilde{y}_y^2', '\tilde{y}_t^5', '\tilde{y}_y^{10}'};

svlabel    = {'pigap', 'ugap', 'intgap'};
for n = 1 : Nyield
    svlabel = cat(2, svlabel, sprintf('yield%dgap', n));
end

datalabel = sprintf('nopinoise.zlb.SPECTREVAR.%s.VARlags2.order1234.ar1svcor.T236.Tjumpoff0', dataset);
filext = sprintf('%s.dat', datalabel);

%% some parameters
fontsize = 14;

yndxPolicyrate = 3;
yndxInflation  = 2;

zlbndx   = yNaN(:,yndxPolicyrate);
startZLB = find(zlbndx, 1, 'first');
stopZLB  = find(zlbndx, 1, 'last');
ELBound  = .25;

fractiles   = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;
ndxmean     = 1;
ndxmedian   = 2;
ndxtails    = 2 + [3 8];
ndxtails2   = 2 + [3 5 6 8];

yearticks   = [datenum(1960:10:2010,1,1) datenum(2018,1,1)];

% END PREAMBLE

