function rndStream = getDefaultStream()
% GETDEFAULTSTREAM ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 23-Dec-2013 12:46:21 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 8.1.0.604 (R2013a) 
% FILENAME  : getDefaultStream.m 


if verLessThan('matlab', '8')
    rndStream = RandStream.getDefaultStream;
else
    rndStream = RandStream.getGlobalStream;
end
