function dspace = loaddebug(loadIntoStruct)
% LOADDEBUG ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 03-Nov-2014 10:34:35 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 8.1.0.604 (R2013a) 
% FILENAME  : loaddebug.m 

if nargin < 1 
    loadIntoStruct = false;
end
if nargout > 0 
    loadIntoStruct = true;
end

% if ~loadIntoStruct
%     evalin('base', 'clear')
% end

dlist = dir('*.debug');
vlist = arrayfun(@(x) strtok(x.name, '.'), dlist, 'UniformOutput', false);
if loadIntoStruct
    for v = 1 : length(vlist)
        dspace.(vlist{v}) = importdata(dlist(v).name);
    end
else
    for v = 1 : length(vlist)
        evalin('base', sprintf('%s = importdata(''%s'');', vlist{v} , dlist(v).name));
    end
end
evalin('base', 'whos')