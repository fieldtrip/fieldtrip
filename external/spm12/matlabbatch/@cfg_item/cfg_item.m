function item = cfg_item(varargin)

% This is the generic configuration item class, from which all other
% classes are derived. 
%
% Data structure
% ==============
% Description fields
%    * name  - display name of config item
%    * tag   - tag of the menu item
%    * val   - (optional) val field: cell array
%    * check - (optional) function handle to implement configuration
%              specific checks based on the harvested subtree rooted at
%              this node. It will be evaluated during harvest if all
%              dependencies in the harvested subtree are resolved and all
%              val's are set. 
%              This function should return an empty string on success and
%              a string explaining why it failed otherwise.
%    * rewrite_job - (optional) function handle to rewrite a job prior to
%              initialisation. This function will be called during
%              initialise(), before any validity checks, subtree
%              initialisation or value assignments are made. The function
%              takes a proposed subjob as input, and should return a valid
%              subjob as output. rewrite_job can be used to implement
%              silent upgrades to jobs when configuration trees have
%              changed.
%    * help  - help text
%    * def   - defaults setting (only evaluated for cfg_leaf items),
%              holding a function handle. This function handle should
%              accept both an empty and a non-empty argument list.
%              If there is no .val{1} set for an cfg_leaf item,
%              feval(def, {}) will be evaluated to retreive a default value.
%              Any value returned that does not match the size/type/filter
%              etc. requirements of the item, will resolve to <UNDEFINED>.
%              To change a default value, feval(def, {newval}) will be
%              called. It is up to the defaults function to decide whether
%              this value will be stored permanently or just for the
%              current instance of the configuration tree. Only values
%              which are valid entries for this field are accepted. If the
%              value is not valid, it will not be changed.
%              To use a registry like defaults function with key/value
%              pairs as arguments, construct the function handle like this:
%              @(defval)get_defaults('some.key', defval{:})
%              This will result in 'get_defaults' being called with the key
%              argument only for retrieving defaults, and with key plus
%              defval arguments to set defaults.
%    * preview - (optional) A function callback that accepts the
%              harvested configuration subtree rooted at this
%              cfg_item. It is evaluated from the GUI and can be used to
%              display information about the entered data. The GUI only
%              calls this callback if the entire subtree is complete
%              (all_leafs/all_set) and contains no dependency objects.
% GUI/job manager fields
%    * expanded
%    * hidden
%
% Public Methods
% ==============
%    * get_strings - returns name of object
%                    No validity check performed here, this needs to be
%                    added in child class method.
%    * gettag      - returns tag
%    * help        - returns help text
%    * harvest     - returns item.val{1}, or '<UNDEFINED>' if empty, see below
%    * all_set     - returns ~isempty(item.val)
%
% Public internal Methods
% =======================
%    * subsasgn
%    * subsref
%    * display
%    * disp
%
% Output in Job Structure (harvest)
% =================================
% cfg_item/harvest returns item.val{1}. If this is a dependency object
% and dependencies shall and can be resolved the contents of the
% dependencies will be returned. Otherwise the dependency object(s) will
% be returned. This is the default behavior for all cfg_leaf items. 
%
% The layout of the configuration tree and the types of configuration items
% have been kept compatible to a configuration system and job manager
% implementation in SPM5 (Statistical Parametric Mapping, Copyright (C)
% 2005 Wellcome Department of Imaging Neuroscience). This code has been
% completely rewritten based on an object oriented model of the
% configuration tree.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_item.m 6133 2014-08-07 10:35:08Z volkmar $

rev = '$Rev: 6133 $'; %#ok

myclass = mfilename;
% Get local fields and defaults from private/mysubs_fields
[fn, defs] = mysubs_fields;
fnd = [fn' defs']';

if nargin == 1
    if isstruct(varargin{1})
        % assume input is a struct to be converted back into a class
        % return with error if this does not work
        if numel(fieldnames(varargin{1})) == numel(fn) && all(isfield(varargin{1}, fn))
            item = class(varargin{1}, myclass);
            return;
        else
            cfg_message('matlabbatch:constructor:reclassify', ['Don''t know how to convert this ' ...
                            'into class ''%s''.'], myclass);
        end;
    end;
    if isa(varargin{1}, myclass)
        item = varargin{1};
        return;
    end;
end;

item = class(struct(fnd{:}), 'cfg_item');
switch nargin
    case 0
        return;
    case 1
        item.name = varargin{1};
    case 2
        item.name = varargin{1};
        item.tag  = varargin{2};
    case 3
        item.name  = varargin{1};
        item.tag   = varargin{2};
        item.check = varargin{3};
    case 4
        item.name  = varargin{1};
        item.tag   = varargin{2};
        item.check = varargin{3};
        item.help  = varargin{4};
    case 5
        item.name         = varargin{1};
        item.tag          = varargin{2};
        item.check        = varargin{3};
        item.help         = varargin{4};
        item.rewrite_job  = varargin{5};
    otherwise
        cfg_message('matlabbatch:constructor:nargin', 'Wrong number of arguments.');
end;
