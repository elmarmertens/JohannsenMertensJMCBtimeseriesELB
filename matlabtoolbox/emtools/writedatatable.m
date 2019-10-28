function writedatatable(wrap, filename, dates, data, datalabels, datefmt)
% WRITEDATATABLE stores times series data in a csv file with column headers
% and a date column.
%
% USAGE: writedatatable(wrap, filename, dates, data, datalabels, datefmt)
% - wrap determines the directory in whic the csv file is to be stored. wrap
%   can be one of the following:
%   -- empty: the current directory (pwd) will be used)
%   -- a string, designating the path to the target directory
%   - a structure with a field called "dir" designating the path to the target directory.
% - filename: name of the target file (without "csv" suffix, which will be appended)
% - dates: Tx1 vector of matlab date numbers
% - data: TxN data matrix (can contain NaNs)
% - datalabels: cell of length N, containing colun labels 
% - datefmt: matlab formatting string for outputting the dates in column 1 
%   (default: % 'yyyyqq')
%  
%   See also: dlmwrite, csvwrite, pwd, datestr

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 28-Mar-2016 10:20:04 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 8.6.0.267246 (R2015b) 
% FILENAME  : writedatatable.m 

if nargin < 6
    datefmt = 'yyyyqq';
end

filename = strcat(filename, '.csv');

if isempty(wrap) 
    return
else
    if isstruct(wrap) && isfield(wrap, 'dir')
        filename = fullfile(wrap.dir, filename);
    else % assuming it is a string ....
        filename = fullfile(wrap, filename);
    end
end

fid = fopen(filename, 'wt');
fprintf(fid, '%15s', 'dates/labels');
fprintf(fid, ',%30s', datalabels{:});
fprintf(fid, '\n');

for n = 1 : length(dates)
    if iscell(dates)
        fprintf(fid, '%15s', dates{n});
    else
        if isempty(datefmt)
            fprintf(fid, '%15d', dates(n));
        else
            fprintf(fid, '%15s', datestr(dates(n), datefmt));
        end
    end
    fprintf(fid, ',%30.16e', data(n,:));
    fprintf(fid, '\n');

    
end
fprintf(fid,'\n');
fclose(fid);

% dlmwrite(filename, data, 'delimiter',',', 'precision', '%30.16e');


     
