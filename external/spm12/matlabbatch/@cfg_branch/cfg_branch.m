function item = cfg_branch(varargin)

% This is the branch configuration item class for non-executable
% branches. It implements branch harvest, all_set, get_strings.
%
% Data structure
% ==============
% Description fields
%    * name  - display name of config item
%    * tag   - tag of the menu item
%    * val   - 1xn cell array of cfg_item objects
%    * check - (optional) function handle to implement configuration
%              specific subsasgn checks based on the harvested subtree
%              rooted at this node
%    * help  - help text
% GUI/job manager fields
%    * expanded
%    * hidden
% All fields are inherited from the generic configuration item class.
%
% Public Methods
% ==============
%    * get_strings - returns name of object
%    * gettag      - returns tag
%    * help        - returns help text
%    * harvest     - returns struct, field names correspond to tags of
%                    items in .val field
%    * all_set     - returns all(all_set(item.val{...}))
%
% Output in Job Structure (harvest)
% =================================
% The resulting structure is a struct. Its fieldnames correspond to the
% tags of the cfg_items in item.val, the value of each field is the
% harvested data of the corresponding child item.
%
% The layout of the configuration tree and the types of configuration items
% have been kept compatible to a configuration system and job manager
% implementation in SPM5 (Statistical Parametric Mapping, Copyright (C)
% 2005 Wellcome Department of Imaging Neuroscience). This code has been
% completely rewritten based on an object oriented model of the
% configuration tree.
%
%   The resulting data structure is a struct, with fieldnames according
%   to the 'tag's of the child nodes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_branch.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

myclass = mfilename;
% Get local fields and defaults from private/mysubs_fields
[fn, defs] = mysubs_fields;

if nargin == 1
    if isstruct(varargin{1})
        % assume input is a struct to be converted back into a class
        % return with error if this does not work
        if numel(fieldnames(varargin{1})) == numel(fn)+2 && ...
                all(isfield(varargin{1}, [fn(:)', {'cfg_item' 'cfg_intree'}]))
            gitem = varargin{1}.cfg_item;
            sitem = rmfield(varargin{1},{'cfg_item','cfg_intree'});
            item  = class(sitem, myclass, gitem, cfg_intree);
            return;
        else
            cfg_message('matlabbatch:constructor:reclassify', ['Don''t know how to convert this ' ...
                                'into class ''%s''.'], myclass);
        end;
    end;
    if isa(varargin{1},myclass)
        item = varargin{1};
        return;
    end;
end;

mxpnargin = 4; % Max 4 arguments to parent initialisation
pnargin = min([nargin,mxpnargin]);
switch nargin
    case 0
        gitem = cfg_item;
    case {1,2,3,4,5}
        gitem = cfg_item(varargin{1:pnargin});
    otherwise
        cfg_message('matlabbatch:constructor:nargin', 'Wrong number of arguments.');
end;
item = class(struct([]), myclass, gitem, cfg_intree);
if nargin > mxpnargin
    item.cfg_item.val{1} = varargin{mxpnargin+1};
    mxpnargin = mxpnargin+1;
end;
