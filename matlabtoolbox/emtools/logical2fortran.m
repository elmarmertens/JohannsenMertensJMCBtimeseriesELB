function logical2fortran(filename, x)
% LOGICAL2FORTRAN writes a matrix of logical values (integer 0/1)
%  
% Usage: logical2fortran(filename, x)
% 
% See also mat2fortran, vec2fortran

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 24-Jun-2012 17:25:16 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : vec2fortran.m 

if islogical(x)
   x = double(x);
else
   warning('em:msg', 'trying to store non-logical array in %s as if it were logical', filename)
end

dlmwrite(filename, x, 'delimiter','', 'precision', '%2d');
