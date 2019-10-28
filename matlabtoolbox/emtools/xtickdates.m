function xtickdates(dates, varargin)
% XTICKDATES sets dateticks 
%  xtickdates(dates) or xtickdates(dates, 10, 'keepticks') etc
%  in addition to calling datetick, the function also implements a bugfix
%  to ensure properly prionted datelabels by disablign 3GL rendering
% see also datetick

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 06-Dec-2012 10:38:58 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.14.0.739 (R2012a) 
% FILENAME  : xtickdates.m 

xlim(dates([1 end]))
datetick('x', 'keeplimits', varargin{:})
% set(gcf, 'Renderer', 'painters')
