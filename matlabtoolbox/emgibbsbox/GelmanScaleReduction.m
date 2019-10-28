function [ScaleReduction, draws, EffectiveDraws] = GelmanScaleReduction(draws, fn, plotFlag)
% GELMANSCALEREDUCTION compute Gelman scale reduction factor (and generate plot)
% [ScaleReduction, MergeDraws] = GelmanScaleReduction(draws, fn, plotFlag, wrap)
%   ...

%   Coded by  Elmar Mertens, em@elmarmertens.com


%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 27-Aug-2009 11:43:53 $
% $Revision : 2.00 $ "Lean and Mean"
% DEVELOPED : 7.7.0.471 (R2008b)
% FILENAME  : GelmanScaleReduction.m


%% pre-amble
if nargin < 2 || isempty(fn)
    fn    = fieldnames(draws);
end

if nargin < 3
    if nargout == 0
        plotFlag = true;
    else
        plotFlag = false;
    end
end

%% fold draws (new PARFOR storage structure)
% new PARFOR structure is draws(n).FIELD
% whereas old structure has draws.FIELD(:,:,n) the latter is more handy
% for computing the Gelman Statisticsacross streams, hence they will be
% used below

draws    = draws(:);
Nstreams = length(draws);
for f = 1 : length(fn)
    celle = arrayfun(@(x) x.(fn{f})(:), draws, 'uniformoutput', false);
    matsize = [size(draws(1).(fn{f})) Nstreams];
    draws = rmfield(draws, fn{f});
    folded.(fn{f}) = cell2mat(celle);
    clear celle
    folded.(fn{f}) = reshape(folded.(fn{f}), matsize);
    
    % OLD CODE:
    % folded.(fn{f}) = reshape(cell2mat(arrayfun(@(x) x.(fn{f})(:), draws, 'uniformoutput', false)), [size(draws(1).(fn{f})) Nstreams]);
    
end
clear draws
draws = folded;
clear folded

%% compute Gelman Stats
for f = 1 : length(fn)
    
    % in Gelman's notation we have 
    %   m = Nstreams (parallel sequences)
    %   n = Nsims    (draws per sequence)
    
    
    % identify the last dimension of the field's matrix,
    % since it contains the streams
    dimms           = ndims(draws.(fn{f}));
    
    % second-to-last contains sweeps per stream
    % (in principle, this should be the same number across fields)
    Nsims           = size(draws.(fn{f}), dimms - 1);
    
    % ACROSS-STREAM VARIANCE OF WITHIN-STREAM MEANS
    % 1) within-stream means
    this            = sum(draws.(fn{f}), dimms - 1) / Nsims; % notice: could compute psi_dj first, then s2 by hand ...
    
    % 2) demean data (used below to compute within stream variance)
    that            = bsxfun(@minus, draws.(fn{f}), this);
    
    % 3) final result
    this            = bsxfun(@minus, this, sum(this, dimms) / Nstreams);
    % notice: Gelman scales B by Nsims, which is however undone again later
    gelmanB         = sum(this.^2, dimms) / (Nstreams - 1); % variance scaled by m-1
    clear this
    
    % ACROSS-STREAM AVERAGE OF WITHIN-STREAM VARIANCES
    gelmanW         = sum(sum(that.^2, dimms - 1), dimms) / Nstreams / (Nsims - 1); 
    % notice the two divisors: Nstreams for the mean computation (outer sum)
    % and (Nsims - 1) for the variance computation (inner sum)
    clear that 
    
    ScaleReduction.(fn{f})  = sqrt((Nsims - 1) / Nsims  + gelmanB ./ gelmanW);
    
    if nargout > 2
        EffectiveDraws.(fn{f})  = Nstreams * (1 + ...
            (Nsims - 1) / Nsims  .* gelmanW ./ gelmanB);
    end
    
    % OLD CODE
    %     this      = draws.(fn{f}); % grab element of draws (notice: this may be wasting a bit of memory)
    %
    %     dimms     = ndims(this); % last dimension contains streams, second-to-last contains sweeps per stream
    %     Nsims     = size(this, dimms - 1); % corresponds to n
    %
    %     s2               = var(this, 0, dimms - 1); % variance scaled by n-1
    %     gelmanW.(fn{f})  = mean(s2, dimms);
    %
    %     psi_dj           = mean(this, dimms - 1); % notice: could compute psi_dj first, then s2 by hand ...
    %     gelmanB.(fn{f})  = var(psi_dj, 0, dimms) * Nsims; % variance scaled by m-1
    %
    %     varplus.(fn{f})         = (Nsims - 1) / Nsims * gelmanW.(fn{f}) + gelmanB.(fn{f}) / Nsims;
    %
    %     ScaleReduction.(fn{f})  = sqrt(varplus.(fn{f}) ./ gelmanW.(fn{f}));
    
end

%% merge draws (if called for)
if nargout > 1
    for f = 1 : length(fn)
        scissor         = size(draws.(fn{f}));
        draws.(fn{f})   = reshape(draws.(fn{f}), [scissor(1:end-2) prod(scissor(end-1:end))]);
    end
end

%% plot (if called for)
if plotFlag
    plotGelmanScaleReduction(ScaleReduction, fn)
end


%% some older code
% 
%             % grab data into dummy variable
%             this      = draws.(fn{f}); % this uses extra memory, but allows us to skip the call to var, where similar extra memory will be used
%             
%             % identify the last dimension of the field's matrix,
%             % since it contains the streams
%             dimms     = ndims(this);
%             
%             % second-to-last contains sweeps per stream (in principle, this should
%             % be the same number across fields
%             nSims     = size(this, dimms - 1);
%             nStreams  = size(this, dimms);
%             
%             % ACROSS-STREAM VARIANCE OF WITHIN-STREAM MEANS
%             % 1) fill gelmanB with means as prelim values
%             gelmanB         = mean(this, dimms - 1); % notice: could compute psi_dj first, then s2 by hand ...
%             % 2) demean data (needed for within stream variance)
%             this            = bsxfun(@minus, this, gelmanB);
%             % 3) final result
%             gelmanB         = var(gelmanB, 0, dimms) * nSims; % variance scaled by m-1
%             
%             % ACROSS-STREAM AVERAGE OF WITHIN-STREAM VARIANCES
%             this            = sum(this.^2, dimms - 1) / (nSims - 1); % assuming real inouts, I omit the call to abs used by Matlab's var functions
%             gelmanW         = mean(gelmanW, dimms);
%             clear this
%             
%             ScaleReduction.(fn{f})  = sqrt(((nSims - 1)  + gelmanB 
