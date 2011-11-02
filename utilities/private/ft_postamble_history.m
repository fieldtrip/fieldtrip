% FT_POSTAMBLE_HISTORY stores the present cfg structure in the output variable
%
% Use as
%   ft_postamble history outputvar

global ft_default

for tmpindx=1:length(ft_default.postamble)
  eval(sprintf('try, %s.cfg = cfg; end', ft_default.postamble{tmpindx}));
end
clear tmpindx

