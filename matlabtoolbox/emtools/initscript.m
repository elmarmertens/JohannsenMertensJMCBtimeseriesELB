%% init some stuff

%   Coded by  Elmar Mertens, em@elmarmertens.com

% clear workspace (but not debug points, nor the profiling state)
clear variables
clear global
% clear functions

% close all figures
close(findobj(allchild(0),'flat','Visible','on'));

fclose all;
clc

format short

lastwarn('', ''); % reset warning state
warning on em:msg

drawnow

tic

%% empty wrap
wrap        = [];
thismfile   = mfilenamecaller;

