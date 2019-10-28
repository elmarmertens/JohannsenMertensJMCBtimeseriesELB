function [rndStream, algorithm] = initRandStream(streamIndex, Nstreams, algorithm, seed)
% INITRANDSTREAMS returns a Random Stream indexed by the (optional) argument streamIndex
% rndStream = initRandStream(streamIndex) returns  
% rndStreams = RandStream.create('mrg32k3a','NumStreams',Nstreams, 'StreamIndices', streamIndex);
%
% usage: [rndStream, algorithm] = initRandStreams(streamIndex, algorithm, seeds)
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
    streamIndex = 1;
end
if nargin < 2 || isempty(Nstreams)
    Nstreams = 1;
end
if nargin < 3 || isempty(algorithm)
    algorithm = '64';
end
if nargin < 4
    seed = sum(clock * 1000);
end

%% prep parameters



switch lower(algorithm)
    case {'mlfg6331_64', '64'}
        algorithm = 'mlfg6331_64';
case {'mrg32k3a', '32'}
        algorithm = 'mrg32k3a';
    otherwise
        error('algorithm %s not recognized', algorithm)
end

%% create RandomStreams
rndStream = RandStream.create(algorithm,'NumStreams', Nstreams, 'StreamIndices', streamIndex, 'Seed', seed);

