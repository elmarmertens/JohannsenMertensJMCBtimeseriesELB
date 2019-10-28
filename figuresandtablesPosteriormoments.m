%% tabualte and plot posterior moments
path(pathdef) % sets path to factory default
addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/

initscript
initwrap

doYGAP = false;

preamble4figurescripts

fontsize = 16;

%% read results

maxlambda       = importdata(fullfile(datadir, sprintf('MAXLAMBDA.%s', filext)));
f               = importdata(fullfile(datadir, sprintf('F.%s', filext)));
shockslopes     = importdata(fullfile(datadir, sprintf('SHOCKSLOPES.%s', filext)));
rbarvar         = loaddat(fullfile(datadir, sprintf('RBARVAR.%s', filext)));
% hvar            = importdata(fullfile(datadir, sprintf('HVAR.%s', filext)));

hrho            = loaddat(fullfile(datadir, sprintf('HRHO.%s', filext)));
hbar            = loaddat(fullfile(datadir, sprintf('HBAR.%s', filext)));

horthrho        = loaddat(fullfile(datadir, sprintf('HORTHRHO.%s', filext)));
horthbar        = loaddat(fullfile(datadir, sprintf('HORTHBAR.%s', filext)));
horthvar        = loaddat(fullfile(datadir, sprintf('HORTHVAR.%s', filext)));

avgtermpremia   = importdata(fullfile(datadir, sprintf('AVGTERMPREMIA.%s', filext)));

% hinno
hSigmavech = loaddat(fullfile(datadir, sprintf('HSIGMA.DRAWS.%s', filext)));

% read priors
prior.rbarvar     = loaddat(fullfile(datadir, sprintf('rbarvar.prior.%s', filext)));
prior.shockslopes = loaddat(fullfile(datadir, sprintf('shockslopes.prior.%s', filext)));
% prior.hvar        = loaddat(fullfile(datadir, sprintf('hvar.prior.%s', filext)));
prior.hrho        = loaddat(fullfile(datadir, sprintf('hrho.prior.%s', filext)));
prior.maxlambda   = loaddat(fullfile(datadir, sprintf('maxlambda.prior.%s', filext)));



%% gap persistence
[~, hpost, hprior] = plotpriorposteriordraws(maxlambda, prior.maxlambda, 0: 0.001 : 1, 'k');
legend([hpost, hprior], 'Posterior', 'Prior', 'location', 'best')
set(gca, 'fontsize', fontsize);
wrapcf(sprintf('%smaxlambda', modellabel), wrap)

%% rbarVOL
% [~, hpost, hprior] = plotpriorposteriordraws(sqrt(rbarvar), sqrt(prior.rbarvar), 0 : .001 : .5, 'k');
% xlim([0 max(xlim)])
% legend([hpost, hprior], 'Posterior', 'Prior', 'location', 'best')
% set(gca, 'fontsize', fontsize);
% wrapcf(sprintf('%srbarvol', modellabel), wrap)

%% rbarVAR
% note: imposes the true analytical prior, not the simulated value
[~, hpost, hprior] = plotpriorposteriordraws(rbarvar, prior.rbarvar, 0 : .001 : .1, 'k', [], .2^2);
xlim([0 max(xlim)])
legend([hpost, hprior], 'Posterior', 'Prior', 'location', 'best')
set(gca, 'fontsize', fontsize);
wrapcf(sprintf('%srbarvar', modellabel), wrap)

%% AVG TP
Ndraws = size(avgtermpremia,2);
% recall: levels of each trend have prior variance of 100 and are uncorrelated; differences thus have variance 200
avgtermpremiaPrior = sqrt(200) * randn(Nyield, Ndraws);
avgtermpremiaPrior = bsxfun(@plus, [2.5;3;3.5], avgtermpremiaPrior);
for i = 1 : Nyield
     [~,hpost,hprior] = plotpriorposteriordraws(avgtermpremia(i,:), avgtermpremiaPrior(i,:), 0 : .0001 : 1.6, 'k');
     legend([hpost, hprior], 'Posterior', 'Prior', 'location', 'best')
     set(gca, 'fontsize', fontsize);
     wrapcf(sprintf('%sTermPremium%d', modellabel, i), wrap)
end


%% B - shockslopes
% convert posterior into B
[Ndraws, Nshockslopes] = size(shockslopes);
B = NaN(Ngap,Ngap,Ndraws);
for n = 1 : Ngap
    B(n,n,:) = 1;
end
offset = 0;
for n = 2 : Ngap
    B(n,1:n-1,:) = permute(shockslopes(:,offset + (1:n-1)), [3 2 1]);
    offset = offset + n - 1;
end

% convert prior into B
[Ndraws, ~] = size(prior.shockslopes);
Bprior = NaN(Ngap,Ngap,Ndraws);
for n = 1 : Ngap
    Bprior(n,n,:) = 1;
end
offset = 0;
for n = 2 : Ngap
    Bprior(n,1:n-1,:) = permute(prior.shockslopes(:,offset + (1:n-1)), [3 2 1]);
    offset = offset + n - 1;
end

midB   = mean(B,3);
lowerB = prctile(B, 5, 3);
upperB = prctile(B, 95, 3);

filename = sprintf('SHOCKSLOPES%s.tex', modellabel);

fid = fopen(filename, 'wt');
% header
fprintf(fid, '\\begin{normalsize}\n');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('c', 1, Ngap)');
fprintf(fid, '\\toprule\n');

fprintf(fid, ' & \\multicolumn{%d}{c}{Choleski Residuals of Shocks to \\ldots} ', Ngap); 
fprintf(fid, '\\\\\n'); % new line
fprintf(fid, '\\cmidrule(lr){%d-%d}\n', 2, 1+Ngap);

fprintf(fid, 'Variables  '); 
fprintf(fid, '& $%s$ ', ynames{:});
fprintf(fid, '\\\\\n'); % new line
fprintf(fid, '\\midrule\n');

for n = 1 : Ngap
    % MID
    fprintf(fid, '$%s$ ', ynames{n});
    for m = 1 : n - 1
        fprintf(fid, '&  {\\bf %6.3f} ', midB(n,m));
    end
    fprintf(fid, '&  {\\bf 1.0}');
    fprintf(fid, '\\\\\n'); % new line
    % TAILS
    for m = 1 : n - 1
        fprintf(fid, '&  [ %6.3f --- %6.3f ] ', lowerB(n,m), upperB(n,m));
    end
    if n == Ngap
        fprintf(fid, '&  \\hspace{2.5cm} '); % todo: this is a little hamfisted
    end
    fprintf(fid, '\\\\\n'); % new line
end    

% footer
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\\end{normalsize}\n');
fclose(fid);
type(filename)

if ~isempty(wrap)
    movefile(filename, wrap.dir)
end
latexwrapper(wrap, 'add', 'sidetab', filename);



%% Construct Correlation-VCV of hSigma
Ndraws    = size(hSigmavech,1);
hSigma    = NaN(Ngap,Ngap,Ndraws);
for n = 1 : Ndraws
     hSigma(:,:,n) = ivechupper(hSigmavech(n,:), [], true);
     hSigma(:,:,n) = hSigma(:,:,n) * hSigma(:,:,n)'; % hSigmavech was constructed from upper-left Choleski factor
     hSigma(:,:,n) = cov2cor(hSigma(:,:,n), true); 
end

Ndraws = 1e5;
hSigmaPrior = iwishdraw(0.2^2 * eye(Ngap), Ngap + 11, Ndraws);
for n = 1 : Ndraws
     hSigmaPrior(:,:,n) = cov2cor(hSigmaPrior(:,:,n), true); 
end

%% tabulate Correlations: Posterior
MIDhSigma   = mean(hSigma,3);
UPPERhSigma = prctile(hSigma,95,3);
LOWERhSigma = prctile(hSigma,5,3);

priorMIDhSigma   = mean(hSigmaPrior,3);
priorUPPERhSigma = prctile(hSigmaPrior,95,3);
priorLOWERhSigma = prctile(hSigmaPrior,5,3);

%% fix simulated priors to identical values 
% prior is symmetric, monte-carlo error introduces unnecessary variations
% just values to first estimate (could also pool across estimates ...)

% diagonals
ndx = eye(Ngap) == 1;
priorMIDhSigma(ndx) = priorMIDhSigma(1,1);
priorUPPERhSigma(ndx) = priorUPPERhSigma(1,1);
priorLOWERhSigma(ndx) = priorLOWERhSigma(1,1);

% off diagonals
priorMIDhSigma(~ndx)   = 0; % priorMIDhSigma(2,1);
taylor = abs(priorUPPERhSigma(2,1));
priorUPPERhSigma(~ndx)  = taylor;
priorLOWERhSigma(~ndx)  = -taylor;

%% start tex file
filename = sprintf('SVcorr%s.tex', modellabel);

fid = fopen(filename, 'wt');
% header
fprintf(fid, '\\begin{normalsize}\n');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('c', 1, Ngap)');
fprintf(fid, '\\toprule\n');

% fprintf(fid, ' & \\multicolumn{%d}{c}{Choleski Residuals of Shocks to \\ldots} ', Ngap); 
% fprintf(fid, '\\\\\n'); % new line
% fprintf(fid, '\\cmidrule(lr){%d-%d}\n', 2, 1+Ngap);

% fprintf(fid, 'Variables  ');
fprintf(fid, '& $%s$ ', ynames{:});
fprintf(fid, '\\\\\n'); % new line
fprintf(fid, '\\midrule\n');

for n = 1 : Ngap
    % MID
    fprintf(fid, '$%s$ ', ynames{n});
    for m = 1 : n - 1
        fprintf(fid, '&  {\\bf %6.3f} ', MIDhSigma(n,m));
    end
    fprintf(fid, '&  {\\bf 1.0}');
    % mid: prior
    for m = n + 1 : Ngap
      fprintf(fid, '&  {\\it %6.3f} ', priorMIDhSigma(n,m));
    end
    fprintf(fid, '\\\\\n'); % new line
    % TAILS
    for m = 1 : n - 1
        fprintf(fid, '&  [ %6.3f --- %6.3f ] ', LOWERhSigma(n,m), UPPERhSigma(n,m));
    end
    fprintf(fid, '&  ');
%     if n == Ngap
%         fprintf(fid, '&  \\hspace{2.5cm} '); % todo: this is a little hamfisted
%     end
    % prior tails
    for m = n + 1 : Ngap
        fprintf(fid, '&  [ {\\it %6.3f --- %6.3f} ] ', priorLOWERhSigma(n,m), priorUPPERhSigma(n,m));
    end
    
    fprintf(fid, '\\\\\n'); % new line
end    

% footer
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
% fprintf(fid,' NOTE: PRIORS ON UPPER HALF ARE SIMULATED, DO SOME ROUNDING OR SO FOR SIMULATION ERROR\n');
fprintf(fid, '\\end{normalsize}\n');
fclose(fid);
type(filename)

if ~isempty(wrap)
    movefile(filename, wrap.dir)
end
if ~isempty(wrap) && wrap.id ~= 0
    latexwrapper(wrap, 'add', 'sidetab', filename);
end

%% Hrho, Hbar, Hvar

hbar = exp(hbar /2);
horthbar = exp(horthbar /2);

% posterior moments
MIDhrho   =  mean(hrho,1);
UPPERhrho =  prctile(hrho,95,1);
LOWERhrho =  prctile(hrho,5,1);

MIDhbar   =  mean(hbar,1);
UPPERhbar =  prctile(hbar,95,1);
LOWERhbar =  prctile(hbar,5,1);

MIDpibarhrho   =  mean(horthrho(:,1));
UPPERpibarhrho =  prctile(horthrho(:,1),95);
LOWERpibarhrho =  prctile(horthrho(:,1),5);

horthvar = sqrt(horthvar);
MIDpibarhvar   =  mean(horthvar(:,1));
UPPERpibarhvar =  prctile(horthvar(:,1),95);
LOWERpibarhvar =  prctile(horthvar(:,1),5);

MIDpibarhbar   =  mean(horthbar(:,1));
UPPERpibarhbar =  prctile(horthbar(:,1),95);
LOWERpibarhbar =  prctile(horthbar(:,1),5);



%% start tex file
filename = sprintf('SVAR1%s.tex', modellabel);
% report rho and sigma of SV AR1's

fid = fopen(filename, 'wt');
% header
fprintf(fid, '\\begin{normalsize}\n');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('c', 1, 3)');
fprintf(fid, '\\toprule\n');

% PANEL A: Trend Inflation
fprintf(fid, '\\multicolumn{%d}{c}{{\\bf PANEL A: Trend Inflation}} ', 4);
fprintf(fid, '\\\\\n'); % new line
fprintf(fid, '\\midrule\n');
fprintf(fid, 'SV in \\ldots \n');
fprintf(fid, ' & \\multicolumn{%d}{c}{$\\exp{(\\mu_{\\bar\\pi} / 2)}$} ', 1); 
fprintf(fid, ' & \\multicolumn{%d}{c}{$\\rho_{\\bar\\pi}$} ', 1); 
fprintf(fid, ' & \\multicolumn{%d}{c}{$\\phi_{\\bar\\pi}$} ', 1); 
fprintf(fid, '\\\\\n'); % new line
fprintf(fid, '\\midrule\n');

% MID
fprintf(fid, '$%s$ ', '\bar\pi_t');
fprintf(fid, '&  {\\bf %6.3f} ', MIDpibarhbar);
fprintf(fid, '&  {\\bf %6.3f} ', MIDpibarhrho);
fprintf(fid, '&  {\\bf %6.3f} ', MIDpibarhvar);
fprintf(fid, '\\\\\n'); % new line
% TAILS
fprintf(fid, '&  [ %6.3f --- %6.3f ] ', LOWERpibarhbar, UPPERpibarhbar);
fprintf(fid, '&  [ %6.3f --- %6.3f ] ', LOWERpibarhrho, UPPERpibarhrho);
fprintf(fid, '&  [ %6.3f --- %6.3f ] ', LOWERpibarhvar, UPPERpibarhvar);
fprintf(fid, '\\\\\n'); % new line
fprintf(fid, '\\midrule\n');

fprintf(fid, '\\multicolumn{%d}{c}{{\\bf PANEL B: Gap Variables}} ', 4);
fprintf(fid, '\\\\\n'); % new line
fprintf(fid, '\\midrule\n');

fprintf(fid, 'SV in \\ldots  ');
fprintf(fid, ' & \\multicolumn{%d}{c}{$\\exp{(\\tilde\\mu_i / 2)}$} ', 1); 
fprintf(fid, ' & \\multicolumn{%d}{c}{$\\tilde\\rho_i$} ', 1); 
fprintf(fid, ' & \\multicolumn{%d}{c}{$\\sqrt{\\tilde\\Phi_{ii}}$} ', 1); 
fprintf(fid, '\\\\\n'); 

fprintf(fid, '\\midrule\n');

for n = 1 : Ngap  
    % MID
    fprintf(fid, '$%s$ ', ynames{n});
    fprintf(fid, '&  {\\bf %6.3f} ', MIDhbar(n));
    fprintf(fid, '&  {\\bf %6.3f} ', MIDhrho(n));
    fprintf(fid, '&  {\\bf %6.3f} ', MIDhSigma(n,n));
    fprintf(fid, '\\\\\n'); % new line
    % TAILS
    fprintf(fid, '&  [ %6.3f --- %6.3f ] ', LOWERhbar(n), UPPERhbar(n));
    fprintf(fid, '&  [ %6.3f --- %6.3f ] ', LOWERhrho(n), UPPERhrho(n));
    fprintf(fid, '&  [ %6.3f --- %6.3f ] ', LOWERhSigma(n,n), UPPERhSigma(n,n));
    fprintf(fid, '\\\\\n'); % new line
end    

% footer
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\\end{normalsize}\n');
fclose(fid);
type(filename)

if ~isempty(wrap)
    movefile(filename, wrap.dir)
end
if ~isempty(wrap) && wrap.id ~= 0
    latexwrapper(wrap, 'add', 'tab', filename);
end


%% finish
finishwrap
dockAllFigures
finishscript
