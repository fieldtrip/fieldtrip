function [script, details] = analysisprotocol(cfg, datacfg)

% ANALYSISPROTOCOL tries to reconstruct the complete analysis protocol that
% was used to create an arbitrary FieldTrip data structure. It will create
% a Matlab script (as text file) and a flowchart with a graphical
% representation.
%
% Use as
%   ANALYSISPROTOCOL(cfg, data)
% where the first cfg input contains the settings that apply to the
% behaviour of this particular function and the second data input argument
% can be the output of any FieldTrip function, e.g. PREPROCESSING,
% TIMELOCKANALYSIS, SOURCEANALYSIS, FREQSTATISTICS or whatever you like.
%
% Alternatively, for the second input argument you can also only give the
% configuration of the processed data (i.e. "data.cfg") instead of the full
% data.
%
% The configuration options that apply to the behaviour of this function are
%  cfg.feedback   = 'no', 'text', 'gui' or 'yes', whether text and/or
%                   graphical feedback should be presented (default = 'yes')
%  cfg.filename   = string, filename of m-file to which the script will be
%                   written (default = [])
%  cfg.remove     = cell-array with strings, determines which objects will
%                   be removed from the configuration prior to writing it to
%                   file. For readibility of the script, you may want to
%                   remove the large objectssuch as event structure, trial
%                   definition, source positions
% cfg.keepremoved = 'yes' or 'no', determines whether removed fields are
%                   completely removed, or only replaced by a short textual
%                   description (default = 'no')
%
% This function uses the nested cfg and cfg.previous that are present in
% the data structure. It will use the configuration and the nested previous
% configurations to climb all the way back into the tree. This funtction
% will print a complete Matlab script to screen (and optionally to file).
% Furthermore, it will show an interactive graphical flowchart
% representation of the steps taken during the analysis. In the flowchart
% you can click on one of the steps to see the configuration details of
% that step.
%
% Note that the nested cfg and cfg.previous in your data might not contain
% all details that are required to reconstruct a complete and valid
% analysis script.

% TODO the output of this function can perhaps be used as input for the wizard function

% Copyright (C) 2006-2009, Robert Oostenveld
%
% $Log: analysisprotocol.m,v $
% Revision 1.4  2009/05/07 07:58:35  roboos
% some general cleanup, don't know the precise details
%
% Revision 1.3  2009/03/23 21:17:53  roboos
% made font slightly larger
%
% Revision 1.2  2009/03/05 09:40:32  roboos
% small change in text feedback
%
% Revision 1.1  2009/03/05 09:38:31  roboos
% renamed cfg2script into analysisprotocol
% created wrapper for deprecated function that gives a warning
% updated documentation
%
% Revision 1.4  2009/03/05 09:07:40  roboos
% added graphical display of the analysis protocol
% changed some details of the recursion, now works with a few persistent variables
%
% Revision 1.3  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.2  2007/05/15 07:02:16  roboos
% updated help
%
% Revision 1.1  2006/12/04 16:03:18  roboos
% new implementation
%

persistent depth   % this corresponds to the vertical   direction in the figure
persistent branch  % this corresponds to the horizontal direction in the figure
persistent parent
persistent info

fieldtripdefs

% set the defaults
if ~isfield(cfg, 'filename'),    cfg.filename    = [];   end
if ~isfield(cfg, 'keepremoved'), cfg.keepremoved = 'no'; end
if ~isfield(cfg, 'feedback'),    cfg.feedback = 'yes';   end

if ~isfield(cfg, 'remove')
  % this is the default list of configuration elements to be removed. These
  % elements would be very large to print and make the script difficult to
  % read. To get a correctly behaving script, you may have to change this.
  cfg.remove = {
    'sgncmb'
    'channelcmb'
    'event'
    'trl'
    'trlold'
    'artfctdef.eog.trl'
    'artfctdef.jump.trl'
    'artfctdef.muscle.trl'
    'pos'
    'inside'
    'outside'
    'grid.pos'
    'grid.inside'
    'grid.outside'
    'version.name'
    'version.id'
    'vol.bnd.pnt'
    'vol.bnd.tri'
    };
end

feedbackgui  = strcmp(cfg.feedback, 'gui') || strcmp(cfg.feedback, 'yes');
feedbacktext = strcmp(cfg.feedback, 'text') || strcmp(cfg.feedback, 'yes');

% we are only interested in the cfg-part of the data
if isfield(datacfg, 'cfg')
  datacfg = datacfg.cfg;
end

% set up the persistent variables
if isempty(depth),  depth = 1; end
if isempty(branch), branch = 1; end

% start with an empty script
script = '';

% get the function call details before they are removed
try
  thisname = getsubfield(datacfg, 'version.name');
  if isstruct(thisname)
    % I think that this was needed for Matlab 6.5
    thisname = thisname.name;
  end
  [p, f] = fileparts(thisname);
  thisname = f;
catch
  thisname = 'unknown';
end

try
  thisid = getsubfield(datacfg, 'version.id');
catch
  thisid = 'unknown';
end

if feedbacktext
  % give some feedback on screen
  fprintf('\n');
  fprintf('recursion depth = %d, branch = %d\n', depth, branch);
  disp(thisname)
  disp(thisid)
end

% remove the fields that are too large or not interesting
for i=1:length(cfg.remove)
  if issubfield(datacfg, cfg.remove{i})
    if feedbacktext
      fprintf('removing %s\n', cfg.remove{i});
    end
    siz = size(getsubfield(datacfg, cfg.remove{i}));
    if strcmp(cfg.keepremoved, 'yes')
      % keep the field, but replace the value with a descriptive string
      datacfg = setsubfield(datacfg, cfg.remove{i}, sprintf('empty - this was cleared by analysisprotocol, original size = [%s]', num2str(siz)));
    else
      datacfg = rmsubfield(datacfg, cfg.remove{i});
    end
  end
end

% convert this part of the configuration to a matlab script
if isfield(datacfg, 'previous')
  thiscfg = rmfield(datacfg, 'previous');
else
  thiscfg = datacfg;
end

code        = printstruct('cfg', thiscfg);
nl          = sprintf('\n');
emptycfg    = sprintf('cfg = [];\n');
sep         = sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
head1       = sprintf('%% version name = %s\n', thisname);
head2       = sprintf('%% version id   = %s\n', thisid);
thisscript  = [nl sep head1 head2 sep emptycfg code nl];

% remember this part of the configuration info
info(branch,depth).name     = thisname;
info(branch,depth).id       = thisid;
info(branch,depth).cfg      = thiscfg;
info(branch,depth).script   = thisscript;
info(branch,depth).this     = [branch depth];
info(branch,depth).parent   = parent;
info(branch,depth).children = {}; % this will be determined later

if isfield(datacfg, 'previous')
  prev   = parent;
  parent = [branch depth]; % this will be used in the recursive call
  if isstruct(datacfg.previous)
    % increment the depth counter
    depth = depth + 1;
    % use recursion to parse the previous section of the tree
    analysisprotocol(cfg, datacfg.previous);
  elseif iscell(datacfg.previous)
    for i=1:length(datacfg.previous(:))
      % increment the branch counter
      branch = branch + (i>1);
      % increment the depth counter
      depth = depth + 1;
      % use recursion to parse each previous section of the tree
      analysisprotocol(cfg, datacfg.previous{i});
    end
  end
  % revert to the orignal parent
  parent = prev;
end

if depth==1
  % the recursion has finished, we are again at the top level

  % the parents were determined while climbing up the tree
  % now it is time to descend and determine the children
  for branch=1:size(info,1)
    for depth=1:size(info,2)
      if ~isempty(info(branch, depth).parent)
        parentbranch = info(branch, depth).parent(1);
        parentdepth  = info(branch, depth).parent(2);
        info(parentbranch, parentdepth).children{end+1} = [branch depth];
      end
    end
  end

  % complement the individual scripts with the input and output variables
  % and the actual function call
  for branch=1:size(info,1)
    for depth=1:size(info,2)
      outputvar = sprintf('var_%d_%d', info(branch, depth).this);
      inputvar = '';
      commandline = sprintf('%s = %s(cfg%s);', outputvar, info(branch, depth).name, inputvar);
      disp(commandline);

    end
  end

  if nargout>0
    % return the complete script as output argument
    script = '';
    for branch=1:size(info,1)
      for depth=1:size(info,2)
        script = [info(branch,depth).script script];
      end
    end
  end

  if nargout>1
    % return the information details as output argument
    details = info;
  end

  if feedbackgui
    fig = figure;
    hold on
    % the axis should not change during the contruction of the arrows,
    % otherwise the arrowheads will be distorted
    axis manual;
    for branch=1:size(info,1)
      for depth=1:size(info,2)
        plotinfo(info(branch,depth));
      end
    end
    axis auto
    axis off
    guidata(fig,info);
    set(fig, 'WindowButtonUpFcn', @button);
    set(fig, 'KeyPressFcn', @key);
  end % feedbackgui

  if ~isempty(cfg.filename)
    % write the complete script to file
    fprintf('writing result to file ''%s''\n', cfg.filename);
    fid = fopen(cfg.filename, 'wb');
    fprintf(fid, '%s', script);
    fclose(fid);
  end

  % clear all persistent variables
  depth  = [];
  branch = [];
  info   = [];
  parent = [];
  fig    = [];
else
  % this level of recursion has finished, decrease the depth
  depth = depth - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotinfo(element)
if isempty(element.name)
  return
end
w = 0.6;
h = 0.3;
% create the 4 courners that can be used as patch
x = [-w/2  w/2  w/2 -w/2];
y = [-h/2 -h/2  h/2  h/2];
% close the patch
x = [x x(end)];
y = [y y(end)];
% move the patch to the desired location
x = element.this(1) + x;
y = element.this(2) + y;

p = patch(x', y', 0);
set(p, 'Facecolor', [1 1 0.6])

label{1} = element.name;
% label{2} = num2str(element.this);

l = text(element.this(1), element.this(2), label);
set(l, 'HorizontalAlignment', 'center')
set(l, 'interpreter', 'none')
set(l, 'fontUnits', 'normalized')
set(l, 'fontSize', 0.03)
% set(l, 'fontName', 'courier')

% draw an arrow to connect this box to its parent
if ~isempty(element.parent)
  base = element.this   - [0 h/2];
  tip  = element.parent + [0 h/2];
  % ARROW(Start,Stop,Length,BaseAngle,TipAngle,Width,Page,CrossDir)
  arrow(base, tip, [], [], [], [], [], []);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = key(h, eventdata, handles, varargin)
% this is just a placeholder for future functionality
% at the moment it does not do anything

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = button(h, eventdata, handles, varargin)
pos = get(get(gcbo, 'CurrentAxes'), 'CurrentPoint');
x = pos(1,1);
y = pos(1,2);
info = guidata(h);
dist = zeros(size(info)) + inf;
for i=1:numel(info)
  if ~isempty(info(i).this)
    dist(i) = norm(info(i).this - [x y]);
  end
end
% determine the box that is nearest by the mouse click
[m, indx] = min(dist(:));
% show the information contained in that box
uidisplaytext(info(indx).script, info(indx).name);

