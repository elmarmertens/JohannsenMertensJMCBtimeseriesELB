%% plot IRFs

path(pathdef) % sets path to factory default
addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/

initscript
wrap = [];

%#ok<*UNRCH>

%% parameters

doIRFunitSV         = true;
doIRFactualObserver = false;
doWrap              = true;

doYgap    = false;

irflags       = 48 + 5 * 4;
irfBaseyear   = 1965;
doAPF         = true;
doJPG         = false;

datadir = '../rollingstoneResults2019/datIRF';

p = 2;

%% load data

titlename = sprintf('particleVARlags%d', p);

if doYgap
    datalabel = 'spectreTB3MSGS020510OutputGapHeadline2018Q4';
    titlename = strcat(titlename, 'Ygap');
else
    datalabel = 'spectreTB3MSGS020510Headline2018Q4';
    titlename = strcat(titlename, 'Ugap');
end

orderCode = '123456';
titlename = strcat(titlename, 'IRFchol');
if doWrap
    initwrap
else
    wrap = [];
end

timestamp  = [];

ELBound = 0.25;

if doYgap
    modellabel = 'baselineYgap';
else
    modellabel = 'baselineUgap';
end


%% get parameters
filext  = sprintf('particles.VARlags%d.%s.ELB2', p, datalabel);

if doIRFunitSV
    filext = strcat(filext, '.irfUnitSV');
end
if doIRFactualObserver
    filext = strcat(filext, '.irfActualObs');
end
filext  = strcat(filext, sprintf('.order%s.dat', orderCode));


if doAPF
    filext  = strcat('APF.', filext);
end



%% get data

y     = importdata(fullfile(datadir, sprintf('YDATA.%s', filext)))';
yNaN  = logical(importdata(fullfile(datadir, sprintf('YNAN.%s', filext))))';
y(yNaN) = NaN;

zlbNdx = yNaN(:,3);
Ylabels = {'cycle', 'pi', 'shadowrate', 'GS02', 'GS05', 'GS10', 'pibar', 'rbar', ...
    'actualrate+', 'actualrate-', 'actual real rate+', 'actual real rate-', ...
    'actualyield2+', 'actualyield2-', ...
    'actualyield5+', 'actualyield5-', ...
    'actualyield10+', 'actualyield10-', ...
    'pigap', ...
    'spread10+', 'spread10-', ...
    'spread10over2+', 'spread10over2-'...
    'actualyieldgap10+', 'actualyieldgap10-', 'yieldgap10'...
    };

Xlabels = {'pibar', 'rbar', 'r2bar', 'r5bar', 'r10bar', ...
    'pigap', 'cycle', 'intgap', 'int2gap', 'int5gap', 'int10gap', ...
    'pigaplag', 'cycle', 'intgaplag', 'int2gaplag','int5gaplag', 'int10gaplag'};
SVlabels = {'pibar', 'pigap', 'cycle', 'intgap', 'int2gap', 'int5gap', 'int10gap'};


yndxCycle       = 1;
yndxPI          = 2;
yndxShadowrate  = 3;
yndxYield10     = 6;
yndxPIbar       = 7;
yndxRbar        = 8;
yndxActualplus  = 9;
yndxActualminus = 10;
yndxRealplus    = 11;
yndxRealminus   = 12;

yndxYield2plus     = 13;
yndxYield2minus    = 14;
yndxYield5plus     = 15;
yndxYield5minus    = 16;
yndxYield10plus    = 17;
yndxYield10minus   = 18;

yndxPIgapminus     = 19;
yndxSpread10plus   = 20;
yndxSpread10minus  = 21;

yndxSpread10over2plus   = 22;
yndxSpread10over2minus  = 23;

yndxyield10gapplus      = 24;
yndxyield10gapminus     = 25;

yndxyield10gap     = 26;

%% define impulses for IRF

ImpulseLabels = {'cycle', 'pi', 'shadowrate', 'yield02', 'yield05'}; % , 'SVEN10'};

Nimpulse = length(ImpulseLabels);

%% read results
T   = size(y,1);
Ny  = size(y,2);
Nyield = 3;
Nbar = 2 + Nyield;
Nx = Nbar + Ny * p;
Nstates = Nbar + Ny; % the states to be shown (counting from the top)
Nsv     = 1 + Ny;


% dates = importdata(fullfile(datadir, sprintf('%s.dates.txt', datalabel)));
dates = importdata(fullfile(pwd, sprintf('%s.dates.txt', datalabel)));
dates = dates(1:T);
datendx = dates >= datenum(irfBaseyear,1,1);

type(fullfile(datadir, strcat('settings.', filext)));

ndxmean     = 1;
ndxmedian   = 2;
ndxtails    = 2 + [3 5 6 8]; % used for the smoother
% fractiles = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;
horizons = [1, 4, 8, 12, 16];

fontsize = 16;

%% LOAD IRF

NNyy = length(Ylabels);
IRFy = NaN(T,irflags+1,NNyy,Nimpulse);
BASELINEy = NaN(T,irflags+1,NNyy);


% load analytical IRF
Nirf = Ny + 2; % analytical irf
for j = 1 : Nirf
    for n = 1 : Nimpulse
        filename = sprintf('IRFy%dimpulse%d.', j, n);
        filename = strcat(filename, filext);
        IRFy(:,:,j,n) = importdata(fullfile(datadir, filename));
    end
end

% add actualrate+
j = Nirf + 1;
for n = 1 : Nimpulse
    filename = sprintf('IRFysim%dimpulse%d.plus.', j, n);
    filename = strcat(filename, filext);
    IRFy(:,:,yndxActualplus,n) = importdata(fullfile(datadir, filename));
end
% add actualrate-
for n = 1 : Nimpulse
    filename = sprintf('IRFysim%dimpulse%d.minus.', j, n);
    filename = strcat(filename, filext);
    IRFy(:,:,yndxActualminus,n) = importdata(fullfile(datadir, filename));
end

% construct real actual rate+ = actual+ - inflation(next)
IRFy(:,1:end-1,yndxRealplus,:) = IRFy(:,1:end-1,yndxActualplus,:) - IRFy(:,2:end,yndxPI,:);

% construct real actual rate- = actual- - (- inflation(next))
IRFy(:,1:end-1,yndxRealminus,:) = IRFy(:,1:end-1,yndxActualminus,:) + IRFy(:,2:end,yndxPI,:);


% add YIELD2+/-
j = Nirf + 2;
for n = 1 : Nimpulse
    filename = sprintf('IRFysim%dimpulse%d.plus.', j, n);
    filename = strcat(filename, filext);
    IRFy(:,:,yndxYield2plus,n) = importdata(fullfile(datadir, filename));
end
for n = 1 : Nimpulse
    filename = sprintf('IRFysim%dimpulse%d.minus.', j, n);
    filename = strcat(filename, filext);
    IRFy(:,:,yndxYield2minus,n) = importdata(fullfile(datadir, filename));
end

% add YIELD2+/-
j = Nirf + 3;
for n = 1 : Nimpulse
    filename = sprintf('IRFysim%dimpulse%d.plus.', j, n);
    filename = strcat(filename, filext);
    IRFy(:,:,yndxYield5plus,n) = importdata(fullfile(datadir, filename));
end
for n = 1 : Nimpulse
    filename = sprintf('IRFysim%dimpulse%d.minus.', j, n);
    filename = strcat(filename, filext);
    IRFy(:,:,yndxYield5minus,n) = importdata(fullfile(datadir, filename));
end

% add YIELD10+/-
j = Nirf + 4;
for n = 1 : Nimpulse
    filename = sprintf('IRFysim%dimpulse%d.plus.', j, n);
    filename = strcat(filename, filext);
    IRFy(:,:,yndxYield10plus,n) = importdata(fullfile(datadir, filename));
end
for n = 1 : Nimpulse
    filename = sprintf('IRFysim%dimpulse%d.minus.', j, n);
    filename = strcat(filename, filext);
    IRFy(:,:,yndxYield10minus,n) = importdata(fullfile(datadir, filename));
end

% construct pigap: pi - pibar
IRFy(:,:,yndxPIgapminus,:) = -(IRFy(:,:,yndxPI,:) - IRFy(:,:,yndxPIbar,:));

% construct spread10+/-
IRFy(:,:,yndxSpread10plus,:)  = IRFy(:,:,yndxYield10plus,:) - IRFy(:,:,yndxActualplus,:);
IRFy(:,:,yndxSpread10minus,:) = IRFy(:,:,yndxYield10minus,:) - IRFy(:,:,yndxActualminus,:);

% construct spread10over2+/-
IRFy(:,:,yndxSpread10over2plus,:)  = IRFy(:,:,yndxYield10plus,:) - IRFy(:,:,yndxYield2plus,:);
IRFy(:,:,yndxSpread10over2minus,:) = IRFy(:,:,yndxYield10minus,:) - IRFy(:,:,yndxYield2minus,:);

% construct yield10gap+/-
IRFy(:,:,yndxyield10gapplus,:)     = IRFy(:,:,yndxYield10plus,:)  - (IRFy(:,:,yndxPIbar,:) + IRFy(:,:,yndxRbar,:));
IRFy(:,:,yndxyield10gapminus,:)    = IRFy(:,:,yndxYield10minus,:) + (IRFy(:,:,yndxPIbar,:) + IRFy(:,:,yndxRbar,:));

IRFy(:,:,yndxyield10gap,:)         = IRFy(:,:,yndxYield10,:) - (IRFy(:,:,yndxPIbar,:) + IRFy(:,:,yndxRbar,:));


IRFy = choppy(IRFy,10);


%% load baseline sims
for j = 1 : yndxYield10minus % no baselines for anything beyond yndxYield10minus
    switch j
        case {yndxActualplus,yndxActualminus}
            thisbase = Ny + 3;
        case {yndxRealplus, yndxRealminus}
            continue
        case {yndxYield2plus, yndxYield2minus}
            thisbase = Ny + 4;
        case {yndxYield5plus, yndxYield5minus}
            thisbase = Ny + 5;
        case {yndxYield10plus, yndxYield10minus}
            thisbase = Ny + 6;
        otherwise
            thisbase = j;
    end
    for n = 1 : Nimpulse
        filename = sprintf('IRFysim%dbaseline.', thisbase);
        filename = strcat(filename, filext);
        BASELINEy(:,:,j) = importdata(fullfile(datadir, filename));
    end
end




%% plot IRF to negative shadow-rate shock for specific points in time
showlags = 24;
% showlags = irflags;
thesedates = datenum([2007; 2009; 2011; 2016],10,1);
linestyle = {'-', '--', '--', '-'};
colors    = {'b', 'r', 'k', 'g'};
hanni = NaN(length(thesedates),1);
n = 3;
for j = [1 2 3 10 14 18 19 23] % [1 2 3, 5, Ny + 4, Ny + 6, Ny + 8, Ny + 12 19 21 23 yndxyield10gap yndxyield10gapminus yndxyield10gapplus]
    figure
    set(gca, 'fontsize', 16)
    hold on
    for t = 1 : length(thesedates)
        thisndx = find(dates == thesedates(t));
        thisIRF = IRFy(thisndx,1:showlags+1,j,n)';
        if j <= Ny
            thisIRF = -thisIRF;
        end
        %         if doYgap && j == 1 % rescale ygap
        %             thisIRF = 4 * thisIRF;
        %         end
        hanni(t) = plot(0:showlags,thisIRF, ...
            'color', colors{t}, 'linewidth', 3, 'linestyle', linestyle{t});
    end
    xlim([0 showlags])
    set(gca, 'xtick', 0:4:showlags)
    
    switch j
        case yndxPI
            ylim([-.3 .125])
        case yndxCycle
            if doYgap
                ylim([-.1 .5])
            else
                ylim([-.1 .5])
            end
        case yndxSpread10over2minus
            if doYgap
                ylim([-.3 0.2])
            else
                ylim([-.3 .2])
            end
        case {14, 18} % longer term yields
            ylim([-.6 0.1])
    end
    
    %     if j == Ny + 6
    %         ylim([-1.2 .2])
    %     end
    plotOrigin
    %     if j == 3
    %         legend(hanni, datestr(thesedates, 10), 'location', 'southeast')
    %         legend('boxoff')
    %     end
    set(gca, 'fontsize', 24)
    wrapcf(sprintf('%sIRFy%dimpulse%dminus',modellabel, j,n), wrap, [], [], [], doJPG)
    %     YLIM = ylim;
    if j ~= 2
        legend(hanni, datestr(thesedates, 10), 'location', 'best')
    else
        legend(hanni, datestr(thesedates, 10), 'location', 'northwest')
    end
    legend('boxoff')
    %     ylim(YLIM)
    wrapcf(sprintf('%sIRFy%dimpulse%dminusWITHLEGEND',modellabel, j,n), wrap, [], [], [], doJPG)
    title(sprintf('Response of %s to %s', upper(Ylabels{j}), upper(ImpulseLabels{n})))
    wrapcf(sprintf('%sIRFy%dimpulse%dminusWITHTITLE',modellabel, j,n), wrap, [], [], [], doJPG)
end

dockAllFigures

%% plot some baselines
showlags = 8;

hanni    = NaN(2,1);

clear baselineset

i = 1;
baselineset(i).label       = 'actual';
baselineset(i).ndxBaseline = yndxActualplus;
baselineset(i).ndxPlus     = yndxActualplus;
baselineset(i).ndxMinus    = yndxActualminus;

i = 2;
baselineset(i).label       = 'yield2';
baselineset(i).ndxBaseline = yndxYield2plus;
baselineset(i).ndxPlus     = yndxYield2plus;
baselineset(i).ndxMinus    = yndxYield2minus;

i = 3;
baselineset(i).label       = 'yield5';
baselineset(i).ndxBaseline = yndxYield5plus;
baselineset(i).ndxPlus     = yndxYield5plus;
baselineset(i).ndxMinus    = yndxYield5minus;

% i = 4;
% baselineset(i).label       = 'yield10';
% baselineset(i).ndxBaseline = yndxYield10plus;
% baselineset(i).ndxPlus     = yndxYield10plus;
% baselineset(i).ndxMinus    = yndxYield10minus;


%% two plot version
for i = 1  % [1 3] % : length(baselineset)
    for t = 224 % [214 224]  % [100 140 180 190 200 220 230]
        for n = 3 % 1 :  Nimpulse
            
            hanni = NaN(3,1);
            
            % actual rate
            % figure
            % figure('pos',[10 10 1800 600])
            
            plusIRF  = IRFy(t,1:showlags+1, baselineset(i).ndxPlus,n)';
            minusIRF = IRFy(t,1:showlags+1, baselineset(i).ndxMinus,n)';
            baseline = BASELINEy(t,1:showlags+1, baselineset(i).ndxBaseline)';
            
            figure
            set(gca, 'fontsize', 18)
            
            hold on
            hanni(1) = plot(0:showlags,baseline, 'k-', 'linewidth', 3);
            hanni(2) = plot(0:showlags,baseline + plusIRF, 'r-', 'linewidth', 3);
            hanni(3) = plot(0:showlags,baseline + minusIRF, 'b-', 'linewidth', 3);
            ylim([0 3])
            xlim([0 showlags])
            % set(gca, 'ytick', [ 0, .25, .5 :.5:10])
            set(gca, 'ytick',  0 :.5:10)
            plothorzline(ELBound, [], 'k--')
            legend(hanni, 'Baseline', 'Positive Impulse', 'Negative Impulse', 'location', 'NorthWest')
            legend('boxoff')
            
            %                 title('FORECAST PATHS')
            
            % orient landscape
            %             wrapcf(sprintf('%sBaselineImpulse%dtime%sLHS', baselineset(i).label, n, datestr(dates(t),'YYYYQQ')), wrap, [], [], [], doJPG)
            wrapcf(sprintf('%s%sBaseline4CASTImpulse%dtime%s', modellabel, baselineset(i).label, n, datestr(dates(t),'YYYYQQ')), wrap, [], [], [], doJPG)
            
            
            figure
            set(gca, 'fontsize', 18)
            
            hold on
            plot(0:showlags,plusIRF, 'r-', 'linewidth', 3);
            plot(0:showlags,minusIRF, 'b-', 'linewidth', 3);
            
            yscale = max(abs(ylim));
            ylim([-yscale yscale])
            xlim([0 showlags])
            plothorzline(0, [], 'k-', 'linewidth', 3)
            
            %                 title('IRF')
            
            % orient landscape
            wrapcf(sprintf('%s%sBaselineIRFImpulse%dtime%s', modellabel, baselineset(i).label, n, datestr(dates(t),'YYYYQQ')), wrap, [], [], [], doJPG)
            
            
        end
    end
end

%% finish
finishwrap
dockAllFigures
finishscript
