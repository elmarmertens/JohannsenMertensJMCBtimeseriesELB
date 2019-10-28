function [rndStreams, algorithm] = initRandStreams(Nstreams, algorithm, seeds)
% INITRANDSTREAMS returns a cell of parallel Random Streams
% rndStreams = initRandStreams(Nstreams) returns a Nstreams x 1 cell  
% where rndStreams{n} = RandStream.create('mrg32k3a','NumStreams',Nstreams, 'StreamIndices', n);
%
% usage: [rndStreams, algorithm] = initRandStreams(Nstreams, algorithm, seeds)
%   ... 

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 27-Aug-2009 10:27:38 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : initRandStreams.m 

%% parse arguments
if nargin < 1
    Nstreams = 1;
end
if nargin < 2 || isempty(algorithm)
    algorithm = '64';
end
if nargin < 3
    seeds = sum(clock * 1000);
end

%% prep parameters
if isscalar(seeds)
    seeds = repmat(seeds, Nstreams, 1);
end
rndStreams  = cell(Nstreams, 1);

switch lower(algorithm)
    case {'mlfg6331_64', '64'}
        algorithm = 'mlfg6331_64';
case {'mrg32k3a', '32'}
        algorithm = 'mrg32k3a';
    otherwise
        error('algorithm %s not recognized', algorithm)
end

%% create RandomStreams
for n = 1 : Nstreams
   rndStreams{n} = RandStream.create(algorithm,'NumStreams', Nstreams, 'StreamIndices', n, 'Seed', seeds(n));
end

