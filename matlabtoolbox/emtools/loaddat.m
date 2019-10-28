function data = loaddat(filename)
% LOADDAT imports a data file or returns empty if file does not exist
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 26-Mar-2013 17:32:11 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.14.0.739 (R2012a) 
% FILENAME  : loaddat.m 

if exist(filename, 'file')
    data = importdata(filename);
    % data = load(filename, 'ascii');
else
    data = [];
end
