% FT_POSTAMBLE_HASTOOLBOX is executed at the end of each FieldTrip
% function to remove other toolboxes that have been added automatically
% by FT_HASTOOLBOX during execution of the specific function.
%
% Use as
%   ft_postamble hastoolbox
%
% See also FT_PREAMBLE, FT_POSTAMBLE, FT_HASTOOLBOX

global ft_default

if ~isempty(ft_default) && isfield(ft_default, 'hastoolbox')
 for i=1:numel(ft_default.hastoolbox)
    ft_info('removing %s from path', ft_default.hastoolbox{i});
    rmpath(genpath(ft_default.hastoolbox{i}));
  end
  ft_default = rmfield(ft_default, 'hastoolbox');
end

