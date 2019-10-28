function savecf(figurename, directory)
%--------------------------------------------------------------
% Prints the current figure to file 'figurename'
%--------------------------------------------------------------
% function savecf(figurename)

%   Coded by  Elmar Mertens, em@elmarmertens.com


if nargin < 2
    directory = [];
end

set(gcf, 'name', figurename)

orient landscape

if ~isempty(directory)
    %            fig = gcf;
    %            fig.PaperUnits = 'inches';
    %            fig.PaperPosition = [.85 9 12.5 9]; % left, bottom width height
    
    %    set(gcf, 'Renderer', 'painters') % to fix date axis bug
    print('-depsc', '-r300', fullfile(directory, figurename));
    
    % print('-dpdf', '-fillpage', '-r300', fullfile(directory, figurename));
    %    print('-dpdf', '-r300', fullfile(directory, figurename));
end