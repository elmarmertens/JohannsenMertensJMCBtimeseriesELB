function poolsize = getparpoolsize()
 

try 
    p = gcp; % get current pool or try to create a new one *if* parallel preferences are set to automatically create a pool when needed.
    % Please see the "Parallel Computing Toolbox" section of the matlab preferences to see whether a pool will be automatically created. 
    % To open these preferences run "preferences('Parallel Computing Toolbox')" from the Matlab command line.
    % If preferences for automatic pool creation are disabled, "gcp" will not attempt to create a new pool
catch poolME
    % one typical cause of problems could be limited availability of licenses
    warning(poolME.identifier, 'There was a problem obtaining a pool of parallel workers; trying to work without one.\n The error message was:\t %s', ...
        poolME.message)
    p = gcp('nocreate'); % If no pool can be created, do not create new one.
end
 
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end