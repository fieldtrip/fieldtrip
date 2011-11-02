% FT_POSTAMBLE_HISTORY stores the present cfg structure in the output variable
%
% Use as
%   ft_postamble history outputvar


global ft_default

% remember the configuration detaild in the output variable
eval(sprintf('try, %s.cfg = cfg; end', ft_default.postamble{1}));
