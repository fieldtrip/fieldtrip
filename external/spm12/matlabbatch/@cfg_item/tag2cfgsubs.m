function [id, stop, rtaglist] = tag2cfgsubs(item, taglist, finalspec, tropts)

% function [id, stop, rtaglist] = tag2cfgsubs(item, taglist, finalspec, tropts)
% Return the index into the values branch of a configuration tree which
% corresponds to a list of tags. 
% This is the generic tag2cfgsubs function, suitable for all leaf
% cfg_items. It stops with success, if the first element in taglist
% matches gettag(item) and item matches finalspec. In this case, it
% returns an empty substruct. If item matches tropts.stopspec or taglist
% has more than one element then stop = true, else stop = false.
% If unsuccessful, it returns an empty cell and stop = true.
% rtaglist contains the remaining tags that were not matched due to a
% stopping criterion.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: tag2cfgsubs.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok
if strcmp(gettag(item), taglist{1}) && match(item, finalspec)
    id = struct('type', {}, 'subs', {});
    stop = numel(taglist) > 1 || (~isempty(tropts.stopspec) ...
                                  && match(item, tropts.stopspec));
else
    id = {};
    stop = true;
end;
rtaglist = taglist(2:end);