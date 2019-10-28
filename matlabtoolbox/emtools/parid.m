function id = parid()
% PARID ...
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 21-Dec-2017 09:25:27 $
% $Revision : 1.00 $
% DEVELOPED : 9.2.0.556344 (R2017a)
% FILENAME  : parid.m

task = getCurrentTask;
if isempty(task)
    id = 1;
else
    id = get(task, 'ID');
end
