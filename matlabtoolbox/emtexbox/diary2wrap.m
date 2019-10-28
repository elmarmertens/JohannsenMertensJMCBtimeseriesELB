function wrap = diary2wrap(diaryname, wrap, includeFlag)
% function wrap = diary2wrap(diaryname, wrap, includeFlag)
%--------------------------------------------------------------
% Starts recording diary listing and makes entry in wrap.
% Note: Need to call "diary off" separately, in order to terminate diary!
%       If necesseary, this happens automatically when calling latexwrapper(wrap, 'compile')
%--------------------------------------------------------------
% or function diary2wrap(wrap, diaryname, includeFlag)
% mandatory is only wrap

%   Coded by  Elmar Mertens, em@elmarmertens.com

% todo: this switches off standard diary w/o further reactivation

narginchk(1,3)

if nargin < 2
   wrap      = diaryname;
   diaryname = [];
elseif isstruct(diaryname) 
   jack      = wrap;
   wrap      = diaryname;
   diaryname = jack;  
end

if ~isstruct(wrap)
   if isempty(wrap)
      return
   else
      error('wrap argument must be struct')
   end
end

if nargin < 3 
    includeFlag = true;
end

if isempty(diaryname)
   [~, jack] = fileparts(tempname);
   diaryname   = sprintf('screen%s.log', jack);
end

if includeFlag
    latexwrapper(wrap, 'add', 'diary', diaryname);
end

if isfield(wrap, 'diary')
   wrap.diary{end+1} = diaryname;
else
   wrap.diary{1} = diaryname;
end

if isfield(wrap, 'dir')
    diaryname = fullfile(wrap.dir, diaryname);
end

diary off

if exist(diaryname, 'file')
    delete(diaryname);
end


diary(diaryname)
disp(datestr(now))

