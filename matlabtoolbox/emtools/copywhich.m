function status = copywhich(fname)
% COPYWHICH copies file from matlab path to current directory
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 27-Aug-2015 15:09:17 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 8.1.0.604 (R2013a) 
% FILENAME  : copywhich.m 


source = which(fname);
s      = copyfile(source, pwd);
if nargout > 0
    status = s;
end
