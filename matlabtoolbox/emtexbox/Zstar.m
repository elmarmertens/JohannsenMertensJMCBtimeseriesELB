function stars = Zstar(Z, onesided)
% function stars = Zstar(Z)

%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 2
    onesided = false;
end

if onesided
    % Note: keeping the abs to allow for answering whether difference is
    % significantly negative or positive
    if  abs(Z) > norminv(0.99, 0, 1)
        stars = '^{\ast\ast\ast}';
    elseif abs(Z) > norminv(0.95, 0, 1)
        stars = '^{\ast\ast}';
    elseif abs(Z) > norminv(0.9, 0, 1)
        stars = '^{\ast}';
    else
        stars = '';
    end
else
    if  abs(Z) > norminv(0.995, 0, 1)
        stars = '^{\ast\ast\ast}';
    elseif abs(Z) > norminv(0.975, 0, 1)
        stars = '^{\ast\ast}';
    elseif abs(Z) > norminv(0.95, 0, 1)
        stars = '^{\ast}';
    else
        stars = '';
    end
end

% same code, using MATLAB stats toolbox
% if  abs(Z) > norminv(0.995) 
%     stars = '^{\ast\ast\ast}';
% elseif abs(Z) > norminv(0.975) 
%     stars = '^{\ast\ast}';
% elseif abs(Z) > norminv(0.95, 0, 1) 
%     stars = '^{\ast}';
% else
%     stars = '';
% end
