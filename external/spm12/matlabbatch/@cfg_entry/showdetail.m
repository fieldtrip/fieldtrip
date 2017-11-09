function str = showdetail(item)

% function str = showdetail(item)
% Display details for a cfg_files item.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: showdetail.m 4867 2012-08-30 13:04:51Z volkmar $

rev = '$Rev: 4867 $'; %#ok

str = showdetail(item.cfg_item);
str{end+1} = 'Class  : cfg_entry';
str{end+1} = 'An input item.';
str = [str; gencode(item.strtype, 'strtype:')];
str = [str; gencode(item.num,     'num    :')];
if ~isempty(item.extras)
    str{end+1} = 'Its contents must match the following additional criteria:';
    str = [str; gencode(item.extras,  'extras :')'];
end