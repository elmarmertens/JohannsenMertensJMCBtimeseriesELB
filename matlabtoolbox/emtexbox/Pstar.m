function stars = Pstar(p)
% function stars = Pstar(Z)

%   Coded by  Elmar Mertens, em@elmarmertens.com

if  p < .01
    stars = '^{\ast\ast\ast}';
elseif p < .05
    stars = '^{\ast\ast}';
elseif p < .1
    stars = '^{\ast}';
else
    stars = '';
end
