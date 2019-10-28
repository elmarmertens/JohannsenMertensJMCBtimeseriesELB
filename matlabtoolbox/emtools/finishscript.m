%% docking
% if usejava('desktop') 
%    dockAllFigures
% end

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% profile running? if yes, report
profilestatus = profile('status');
if usejava('jvm') && strcmp(profilestatus.ProfilerStatus, 'on')
   profile report
   profile off
end

%% close all files
fclose all;
