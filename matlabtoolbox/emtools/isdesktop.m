function id = isdesktop
% ISDESKTOP returns true if the MATLAB desktop is running
%  
%   isdesktop = usejava('desktop')
% 
% see also usejava

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 29-Mar-2011 17:14:39 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.11.0.584 (R2010b) 
% FILENAME  : isdesktop.m 



id = usejava('desktop');
