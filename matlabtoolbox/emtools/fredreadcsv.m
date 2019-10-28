function DATA = fredreadcsv(fredid, label, SourceFrequency, TargetFrequency, Conversion, plotFlag)
% function DATA = fredreadcsv(fredid, label, SourceFrequency, TargetFrequency, Conversion, plotFlag)
% Possible Conversion methods: "first", "average", "last"
% OUTPUT: DATA.VALUES, DATA.DATES (matlab datenumbers), DATA.NAME = label
% new version, supporting new file structure of FRED. Assumes that

%   Coded by  Elmar Mertens, em@elmarmertens.com

% updated from fredread in Jan 2017 to handle csv format (txt is
% deprecated)

% Elmar Mertens
% Study Center Gerzensee
% elmar.mertens@szgerzensee.ch

% narginchk(1, 5)


if nargin < 2
    label = [];
end

if nargin >= 3 && ~isempty(SourceFrequency)
    SourceFrequency = lower(SourceFrequency);
    switch SourceFrequency
        case 'q'
            SourceFrequency = 'quarterly';
        case 'm'
            SourceFrequency = 'monthly';
        case {'y', 'a'}
            SourceFrequency = 'yearly';
    end
else
    SourceFrequency = 'quarterly';
end

if nargin >= 4 && ~isempty(TargetFrequency)
    TargetFrequency = lower(TargetFrequency);
    switch TargetFrequency
        case 'q'
            TargetFrequency = 'quarterly';
        case 'm'
            TargetFrequency = 'monthly';
        case {'y', 'a'}
            TargetFrequency = 'yearly';
    end
else
    TargetFrequency = [];
end

if nargin < 5 || isempty(Conversion)
    Conversion = 'first';
else
    Conversion = lower(Conversion);
end

if nargin < 6
    plotFlag = 0;
end

fredid   = upper(fredid);
fredfile = strcat(fredid, '.csv');
[fid, msg] = fopen(fredfile, 'r');
if fid == -1
    error('error opening %s:\n%s', fredfile, msg);
end

% S = [];
% R = [];

foo       = importdata(sprintf('%s.csv', fredid));
DATES   = datenum(foo.textdata(2:end,1));
VALUES  = foo.data;
% FREDID  = fredid;
% NAME    = fredid;


% while ~strcmp(S, 'DATE') && ~strcmp(R, 'VALUE')
%     PARSE = fgetl(fid);
%     [S, R] = strtok(PARSE); % % minor TODO: read out FRED titles & ID
%     R = strtrim(R);
%     if strcmpi(S, 'Frequency:')
%         SourceFrequency = lower(R);
%     end
% end
% 
% PARSE  = textscan(fid, '%4n%2n%2n %[1234567890.-]', 'delimiter', '-');
% % old call of textscan below (the new call allows for empty values, which can happen with daily data)
% % PARSE  = textscan(fid, '%4n%2n%2n %f', 'delimiter', '-');


% VALUES = NaN * ones(size(PARSE{4}));
% for t = 1 : length(PARSE{4})
%     tmp = str2double(PARSE{4}{t});
%     if ~isempty(tmp)
%         VALUES(t) = tmp;
%     end
% end



% DATES  = datenum([PARSE{1} PARSE{2} PARSE{3}]);
[Y, M, D] = datevec(DATES); %#ok<ASGLU>
% if ~all(Y == PARSE{1}) || ~all(M == PARSE{2}) || ~all(D == PARSE{3})
%     error('something wrong with date I/O')
% end

% TargetFrequency conversion, if necessary
if strcmp(TargetFrequency, SourceFrequency) || isempty(TargetFrequency)
    TargetFrequency = SourceFrequency;
    % nothing to convert
else
    
    switch SourceFrequency
        case 'monthly'
            switch TargetFrequency
                case 'quarterly'
                    Q      = floor((M - 1) / 3) + 1;
                    
                    switch Conversion
                        case 'first'
                            ndx     = (M == 1) | (M == 4) | (M == 7) | (M == 10);
                            VALUES  = VALUES(ndx);
                            DATES   = DATES(ndx);
                        case 'last'
                            ndx     = (M == 3) | (M == 6) | (M == 9) | (M == 12);
                            VALUES  = VALUES(ndx);
                            ndx     = (M == 1) | (M == 4) | (M == 7) | (M == 10);
                            DATES   = DATES(ndx); % use beginning of quarter dates
                        case 'avg'
                            YQ      = Y + Q / 10;
                            % set data for incomplete quarters to NaN
                            % in priciple, should be sufficient to check for initial/final months 
                            % but, to be safe, let's loop over the whole data vector
                            for t = 1 : length(YQ)
                                if sum(YQ == YQ(t)) ~= 3
                                    VALUES(t) = NaN;
                                end
                            end
                            VALUES  = averageValues(VALUES, YQ);
                            ndx     = (M == 1) | (M == 4) | (M == 7) | (M == 10); % for congruence with FRED
                            DATES   = DATES(ndx);
                        otherwise
                            error('required conversion method unknown');
                    end % switch monthly:quarterly:TargetFrequency
                    
                case 'yearly'
                    warning('em:msg', 'monthly to annual freq conversion is beta version')
                    
                    switch Conversion
                        case 'first'
                            ndx     = (M == 1);
                            VALUES  = VALUES(ndx);
                            DATES   = DATES(ndx);
                        case 'last'
                            ndx     = (M == 12);
                            VALUES  = VALUES(ndx);
                            DATES   = DATES(ndx);
                        case 'avg'
                            VALUES  = averageValues(VALUES, Y);
                            ndx     = (M == 1); % this timing is congruent with FRED
                            DATES   = DATES(ndx);
                        otherwise
                            error('required conversion method unknown');
                    end % switch monthly:annual:TargetFrequency
                    
                otherwise
                    error('other target frequencies currently not convertible')
            end % switch monthly:TargetFrequency
        case {'daily', 'daily, seven day'}
            
            switch TargetFrequency
                case 'monthly'
                    
                    switch Conversion
                        case {'first', 'last'}
                            
                            nanny = ~isnan(VALUES);
                            %                      ndx = ~logical(ones(size(DATES)));
                            ndx = false(size(DATES));
                            [nextyear, nextmonth] = datevec(DATES(1));
                            
                            if strcmp(Conversion, 'last')
                                tau = find(Y == nextyear & M == nextmonth & ~isnan(VALUES), 1, 'last');
                            else
                                tau = 1;
                            end
                            ndx(tau) = true;
                            
                            while ~isempty(tau);
                                if nextmonth  < 12
                                    nextmonth = nextmonth + 1;
                                else
                                    nextmonth = 1;
                                    nextyear  = nextyear + 1;
                                end
                                
                                tau = find(Y == nextyear & M == nextmonth & nanny, 1, Conversion);
                                if ~isempty(tau)
                                    ndx(tau) = true;
                                end
                            end
                            
                            lolly = find(ndx, 1, 'last');
                            if M(lolly) ~= M(end) || Y(lolly) ~= Y(end)
                                warning('em:debug', 'something fishy ...')
                                keyboard
                            end
                            
                            VALUES  = VALUES(ndx);
                            
                            % picking the first day of the month matches the convention used by FRED in its monthly time series
                            DATES  = unique(datenum(Y, M, 1));
                            
                        case 'avg'
                            nanny  = ~isnan(VALUES); % FRED's daily tapes may have missing values
                            
                            % picking the first day of the month matches the convention used by FRED in its monthly time series
                            DATES           = datenum(Y(nanny), M(nanny), 1);
                            [VALUES, DATES] = averageValues(VALUES(nanny), DATES);
                            
                            % VALUES = grpstats(VALUES(nanny), {Y(nanny) M(nanny)});
                            
                            
                        otherwise
                            error('required conversion method unknown')
                    end
                    
                case 'quarterly'
                    Q      = floor((M - 1) / 3) + 1;
                    
                    switch Conversion
                        %                         case 'first'
                        %                             ndx     = (M == 1) | (M == 4) | (M == 7) | (M == 10);
                        %                             VALUES  = VALUES(ndx);
                        %                             %                             DATES   = DATES(ndx);
                        %                         case 'last'
                        %                             ndx     = (M == 3) | (M == 6) | (M == 9) | (M == 12);
                        %                             VALUES  = VALUES(ndx);
                        %                             %                             ndx     = (M == 1) | (M == 4) | (M == 7) | (M == 10);
                        %                             %                             DATES   = DATES(ndx); % use beginning of quarter dates
                        case 'avg'
                            YQ      = Y + Q / 10;
                            VALUES  = averageValues(VALUES, YQ);
                            
                        otherwise
                            error('required conversion method unknown');
                    end % switch daily:quarterly:TargetFrequency
                    
                     DATES   = unique(datenum(Y, M, 1));
                     M       = month(DATES);
                     ndx     = (M == 1) | (M == 4) | (M == 7) | (M == 10); % for congruence with FRED
                     DATES   = DATES(ndx);
                     
                otherwise
                    error('other target frequencies currently not convertible')
            end % switch daily:TargetFrequency
            
        case 'quarterly'
            switch TargetFrequency
                case 'monthly'
                    QDATES = DATES;
                    QVALUES = VALUES;
                    switch Conversion
                        case 'last'  % store at last month of the quarter
                            % 1 - create monthly date vector
                            MDATES = genrMdates(year(QDATES(1)), year(QDATES(end)), 1);
                            MDATES = MDATES(MDATES >= QDATES(1) & MDATES <= QDATES(end));
                            % 2 - shift dates to the end of the quarter
                            [Y, M, D] = datevec(MDATES);
                            DATES = datenum(Y,M+2,D); % note: M is always first month of quarter and D is always 1 (FRED convention)
                            % 3 - match dates and values
                            VALUES = NaN(size(DATES));
                            ndx = ismember(MDATES, QDATES); % Notiece: MDATES and DATES have the same number of elements!
                            VALUES(ndx) = QVALUES;
                        case 'first' % store at first month of the quarter
                            % 1 - create monthly date vector
                            MDATES = genrMdates(year(QDATES(1)), year(QDATES(end)), 1);
                            DATES  = MDATES(MDATES >= QDATES(1) & MDATES <= QDATES(end));
                            % 2 - match dates and values
                            VALUES = NaN(size(DATES));
                            ndx = ismember(DATES, QDATES);
                            VALUES(ndx) = QVALUES;
                        otherwise
                            error('required conversion method unknown');
                    end % switch quarterly:monthly:TargetFrequency
                    
                case 'yearly'
                    warning('em:msg', 'quarterly to annual freq conversion is beta version')
                    switch Conversion
                        case 'first'
                            ndx     = (M <= 3);    % not relying on FRED always putting Jan = Q1
                            VALUES  = VALUES(ndx);
                            DATES   = DATES(ndx);
                        case 'last'
                            ndx     = (M >= 10);
                            VALUES  = VALUES(ndx);
                            DATES   = DATES(ndx);
                        case 'avg'
                            VALUES  = averageValues(VALUES, Y);
                            ndx     = (M == 1); % for congruence with FRED
                            DATES   = DATES(ndx);
                        otherwise
                            error('required conversion method unknown');
                    end % switch quarterly:annual:TargetFrequency
                otherwise
                    error('Conversion of quarterly data to "%s" currently not supported', TargetFrequency)
            end % switch quarterly:TargetFrequency
            
        otherwise
            error('Data from source frequency %s currently not convertible', SourceFrequency)
    end % switch SourceFrequency
    
end

if isempty(label)
    label = fredid;
end

% copy the final resuls into output structure
DATA = struct('VALUES', VALUES, 'DATES', DATES, 'NAME', label, 'FREDID', fredid);
% DATA.VALUES = VALUES;
% DATA.DATES  = DATES;
% DATA.NAME   = label;
% DATA.FREDID = fredid;

fclose(fid);


% gimmick: an instant plot, if required
if plotFlag
    figure
    plot(DATA.DATES, DATA.VALUES);
    title(sprintf('%s (%s) [%s->%s %s]', label, fredid, SourceFrequency, TargetFrequency, Conversion))
    switch lower(TargetFrequency)
        case 'daily'
            datetick('x', 25)
        case 'monthly'
            datetick('x', 12)
        case 'quarterly'
            datetick('x', 17)
        otherwise
            error('TargetFrequency not recognized')
    end
end


function [avgVALUES, uniqueTags] = averageValues(VALUES, tag)
% converts monthly values into quarterly averages

uniqueTags = unique(tag);

avgVALUES = NaN(length(uniqueTags), 1);
for n = 1 : length(uniqueTags)
    thisquarter     = tag == uniqueTags(n);
    avgVALUES(n)    = mean(VALUES(thisquarter));
end

