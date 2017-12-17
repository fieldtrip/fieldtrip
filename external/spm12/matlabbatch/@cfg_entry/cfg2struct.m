function sitem = cfg2struct(item)

% function sitem = cfg2struct(item)
% Return a struct containing all fields of item plus a field type. This is
% the method suitable for entry classes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg2struct.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

% Get parent struct, re-classify as field 'type'
sitem = cfg2struct(item.cfg_item);
sitem.type = class(item);

% Need to cycle through added fields
fn = mysubs_fields;
for k = 1:numel(fn)
    sitem.(fn{k}) = item.(fn{k});
end;
