function [item, defaults] = val2def(item, defaults, funname, deftag)
% function [item, defaults] = val2def(item, defaults, funname, deftag)
% If a cfg_leaf item has a value, extract it and generate code for defaults
% retrieval. This function works in a way similar to harvest, but with a
% much simpler logic. Also, it modifies the returned configuration tree by
% clearing the .val fields if they are moved to defaults. If a .def field
% is already present, its value will be evaluated and a new callback will
% be installed.
% Initially, defaults and deftag should be empty.
% This function is identical for all cfg_leaf classes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: val2def.m 3469 2009-10-16 08:43:15Z volkmar $

rev = '$Rev: 3469 $'; %#ok

if isempty(item.def)
    if ~isempty(item.val) && ~isa(item.val{1},'cfg_dep')
        % Read val entry
        defval = item.val{1};
    end
else
    % Read defaults entry
    defval = item.def({});
end
if exist('defval','var')
    % Create defaults entry
    evalc(['defaults.' deftag ' = defval;']);
    % Create callback
    evalc(sprintf('item.def = @(val)%s(''%s'', val{:});', funname, deftag));
end
% Clear val
item.val = {};
