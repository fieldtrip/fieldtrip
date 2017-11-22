function str = showdoc(item, indent)

% function str = showdoc(item, indent)
% Generic showdoc function for cfg_item classes. It displays the
% (indented) name of the item and the justified help text for this item.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: showdoc.m 3944 2010-06-23 08:53:40Z volkmar $

rev = '$Rev: 3944 $'; %#ok

if isempty(indent)
    str{1} = item.name;
else
    str = {sprintf('%s %s', indent, item.name)};
end
if ~isempty(item.help)
    str = [str(:); item.help(:)]';
end
if ~isempty(item.def)
    str{end+1} = sprintf(['This item has a default value, set via a call ' ...
                        'to function']);
    str{end+1} = func2str(item.def);
end
