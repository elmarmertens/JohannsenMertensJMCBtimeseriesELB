%% plot predictive densities for baeline and zeroesasdata at different points in time
% dataset1: LHS plots of baseline
% dataset1: RHS plots when ELB is treated as data (zeroesasdata option)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/

%% identify where results are store

% datalabels used in estimation
dataset1   = 'spectreTB3MSGS020510OutputGapHeadline2018Q4';
dataset2   = 'spectreTB3MSGS020510OutputGapHeadline2018Q4';

% directories with results files (can be identical)
datadir1 = datadirfn();
datadir2 = datadirfn();


%% general parameters
NNyy    = 8;
ELBound = 0.25;

dates = genrQdates(1960,2030); % importdata(fullfile(datadir1, sprintf('%s.dates.txt', dataset1)));

fractiles   = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;
ndxmean     = 1;
ndxmedian   = 2;
ndxtails    = 2 + [3 5 6 8];

ypredLabel    = {'cycle', 'inf', 'i', 'yield2', 'yield5',  'yield6', 'inf-4q', 's'};
fontSize      = 16;
ylimitswide   = [-5 6.5]; % yaxis limits; set to empty to adapt automatically
ylimitsnarrow = [-2 4.5]; % yaxis limits; set to empty to adapt automatically
theseYticks   = [-6 -4 -2 ELBound 2 4 6];
fixedXlim     = false;

Tlist = [196 197 200 204 224 228];

for n = 8 % plot only short rate forecasts (third element in prediction vector)
    
    for t = 1 : length(Tlist)
        
        T = Tlist(t);
        
        if T > 200
            ylimits = ylimitsnarrow;
        else 
            ylimits = ylimitswide;
        end
        
        datalabel1 = sprintf('nopinoise.zlb.SPECTREVAR.%s.VARlags2.order1234.ar1svcor.T%d.Tjumpoff0', dataset1, T);

        filext1 = sprintf('%s.dat', datalabel1);
        
        %% load dataset1
        datadir   = datadir1;
        filext    = filext1;
        
        YHORIZONS = loaddat(fullfile(datadir, sprintf('YHORIZONS.%s', filext)));
        YPRED     = loaddat(fullfile(datadir, sprintf('YPRED.%s', filext)));
        % NNyy       = size(YPRED,1);
        
        YNANPRED  = logical(loaddat(fullfile(datadir, sprintf('YNANPRED.%s', filext))));
        YPRED(YNANPRED) = NaN;
        
        YDRAW1 = NaN(max(YHORIZONS), 12, NNyy);
        for m = 1 : NNyy
            YDRAW1(:,:,m)  = importdata(fullfile(datadir, sprintf('YDRAW%d.%s', n, filext)));
        end
        YDRAW1observed  = importdata(fullfile(datadir, sprintf('YDRAW%d.%s', 3, filext)));
        

        %% plot predictive density -- baseline
        %         hanni = NaN(3,1);
        xdates = dates(T + (1:max(YHORIZONS)));
        figure
        set(gca, 'fontsize', fontSize)
        hold on
        plotCIybase(YDRAW1(:,ndxmean,n), YDRAW1(:,ndxtails,n), xdates, ylimits(1), 'k-', 'linewidth', 4);
        
        % fill in shadow rate
        a1 = plot(xdates, YDRAW1(:,ndxmean), 'k-', 'linewidth', 4);
        a2 = plot(xdates, YDRAW1observed(:,ndxmean), 'y-.', 'linewidth', 4);
        plotOrigin('k-', [], ELBound)
        if t==1
            legend([a1,a2], 'Shadow rate (mean)', 'Actual rate (mean)', 'Location', 'NorthWest');
            legend boxoff
        end
        if fixedXlim
            theseDates  = dates([193 end]); %#ok<*UNRCH>
            theseXticks = dates(193 : 4 : end); 
        else
            theseDates  = dates(T + YHORIZONS([1 end]));
            theseXticks = dates(T + YHORIZONS([1 3 5 7])); 
        end
        xlim(theseDates)
        if ~isempty(ylimits)
            ylim(ylimits)
        end
        xticks(theseXticks)
        yticks(theseYticks)
        datetick('x', 27, 'keeplimits', 'keepticks')
        set(gca, 'box', 'off')
        set(gcf, 'Renderer', 'painters') % to fix date axis bug     
        print('-djpeg', '-r500', ['SHADOWRATEdensityT',num2str(Tlist(t))]);
    end
end
