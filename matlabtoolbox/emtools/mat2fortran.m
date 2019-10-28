function mat2fortran(filename, x)
% MAT2FORTRAN writes a matrix in a text file for import into fortran 
%  
% Usage: mat2fortran(filename, x)
% 
% See also vec2fortran, logical2fortran

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 24-Jun-2012 17:25:16 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : vec2fortran.m 

dlmwrite(filename, x, 'delimiter','', 'precision', '%30.16e');