function wrapcf(figurename, wrap, captionname, figurecomment, landscape, doJPG)
%--------------------------------------------------------------
% Prints the current figure to file 'figurename' as fig, eps, jpg and pdf.
% inserts into wrap
%--------------------------------------------------------------
% function wrapcf(figurename, wrap, captionname, figurecomment, landscape)

%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 3
   captionname = [];
end
if nargin < 4
   figurecomment = [];
end
if nargin < 5
   landscape = false;
end
if nargin < 6
   doJPG = false;
end

if ~isempty(figurename)
    % replace any "." by "-" in figurename, otherwise, print will not append
    % the appropriate file name extension
    figurename = strrep(figurename, '.', '-');
end

if ~isempty(captionname)
   set(gcf, 'name', captionname)
elseif ~isempty(figurename)
   set(gcf, 'name', figurename)
end

if isempty(wrap) % do nothing, except changing the name of the figure as above
    return
end

% set(gcf, 'Renderer', 'painters') % to fix date axis bug

if nargin > 1 && ~isempty(wrap)
   if (wrap.id ~= 0)
      if landscape
         latexwrapper(wrap, 'add', 'figure', figurename, captionname, figurecomment);
      else
         latexwrapper(wrap, 'add', 'sidewaysfigure', figurename, captionname, figurecomment);
      end
   else
      warning('em:msg', 'Cannot append figure to latexwrapper since wrap file is closed (fid=0)')
   end
   if isfield(wrap, 'dir')
      figurename = fullfile(wrap.dir, figurename);
   end
end


if landscape
   orient landscape
end

if doJPG
    print('-djpeg', '-r500', figurename);
else
        print('-depsc', '-r300', '-loose', figurename);
    %     orient landscape
    %     print('-dpdf', '-r300', '-fillpage', figurename);
end
