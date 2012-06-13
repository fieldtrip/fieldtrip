function [script, details] = ft_analysisprotocol(cfg, datacfg)

% FT_ANALYSISPROTOCOL tries to reconstruct the complete analysis protocol that
% was used to create an arbitrary FieldTrip data structure. It will create
% a Matlab script (as text file) and a flowchart with a graphical
% representation.
%
% Use as
%   ft_analysisprotocol(cfg, data)
%
% where the first cfg input contains the settings that apply to the
% behaviour of this particular function and the second data input argument
% can be the output of any FieldTrip function, e.g. FT_PREPROCESSING,
% FT_TIMELOCKANALYSIS, FT_SOURCEANALYSIS, FT_FREQSTATISTICS or whatever you like.
%
% Alternatively, for the second input argument you can also only give the
% configuration of the processed data (i.e. "data.cfg") instead of the full
% data.
%
% The configuration options that apply to the behaviour of this function are
%  cfg.feedback   = 'no', 'text', 'gui' or 'yes', whether text and/or
%                   graphical feedback should be presented (default = 'yes')
%  cfg.showinfo   = string or cell array of strings, information to display
%                   in the gui boxes, can be any combination of
%                   'functionname', 'revision', 'matlabversion',
%                   'computername', 'username', 'calltime', 'timeused',
%                   'memused', 'workingdir', 'scriptpath' (default =
%                   'functionname', only display function name). Can also
%                   be 'all', show all info. Please note that if you want
%                   to show a lot of information, this will require a lot
%                   of screen real estate.
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
%
% See also FT_PREPROCESSING,FT_TIMELOCKANALYSIS, FT_SOURCEANALYSIS, FT_FREQSTATISTICS

% TODO the output of this function can perhaps be used as input for the wizard function

% Copyright (C) 2006-2009, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

persistent depth   % this corresponds to the vertical   direction in the figure
persistent branch  % this corresponds to the horizontal direction in the figure
persistent parent
persistent info

revision = '$Id$';

% callinfo feedback is highly annoying in this recursive function
% do this here, otherwise ft_defaults will override our setting
if ~isfield(cfg, 'showcallinfo'), cfg.showcallinfo = 'no';   end

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% set the defaults
if ~isfield(cfg, 'filename'),    cfg.filename    = [];   end
if ~isfield(cfg, 'showinfo'),    cfg.showinfo    = {'functionname'};   end
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
    'vol.bnd.pnt'
    'vol.bnd.tri'
    };
elseif ~iscell(cfg.remove)
  cfg.remove = {cfg.remove};
end

if strcmp(cfg.showinfo, 'all')
  cfg.showinfo = {
    'functionname', 'revision', 'matlabversion',...
    'computername', 'username', 'calltime', 'timeused',...
    'memused', 'workingdir', 'scriptpath'...
  };
end

if ~isfield(cfg, 'showinfo')
  cfg.showinfo = {'functionname'};
elseif ~iscell(cfg.showinfo)
  cfg.showinfo = {cfg.showinfo};
end

feedbackgui  = strcmp(cfg.feedback, 'gui') || strcmp(cfg.feedback, 'yes');
feedbacktext = strcmp(cfg.feedback, 'text') || strcmp(cfg.feedback, 'yes') || strcmp(cfg.feedback, 'verbose');
feedbackverbose = strcmp(cfg.feedback, 'verbose');

% we are only interested in the cfg-part of the data
if isfield(datacfg, 'cfg')
  datacfg = datacfg.cfg;
end

% set up the persistent variables
if isempty(depth),  depth = 1; end
if isempty(branch), branch = 1; end

if depth==1 && branch==1 && feedbacktext
  fprintf('steps found in analysis pipeline:\n\n');
end

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

if feedbackverbose
  % give some feedback on screen
  fprintf('\n');
  fprintf('recursion depth = %d, branch = %d\n', depth, branch);
  disp(thisname)
  disp(thisid)
elseif feedbacktext
  % give abridged feedback
  fprintf('%-30s  depth = %2d  branch = %2d\n', thisname, depth, branch);
end

% remove the fields that are too large or not interesting
for i=1:length(cfg.remove)
  if issubfield(datacfg, cfg.remove{i})
    if feedbackverbose
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

% this will keep track of whether an extra branch has been found in this
% step
extraBranchesFound = 0;

prev   = parent;
parent = [branch depth]; % this will be used in the recursive call
  
if isfield(datacfg, 'previous')
  
  if isstruct(datacfg.previous)
    % single previous cfg, no branching here
    
    % increment the depth counter
    depth = depth + 1;
    
    % increment the branch counter
    branch = branch + extraBranchesFound;
    extraBranchesFound = 1;
    
    % use recursion to parse the previous section of the tree
    ft_analysisprotocol(cfg, datacfg.previous);
    
  elseif iscell(datacfg.previous)
    % multiple previous cfgs, branch
    
    for i=1:length(datacfg.previous(:))
      % increment the depth counter
      depth = depth + 1;

      % increment the branch counter
      branch = branch + extraBranchesFound;
      extraBranchesFound = 1;
      
      % use recursion to parse each previous section of the tree
      ft_analysisprotocol(cfg, datacfg.previous{i});
      
    end
  end
end

% check all fields of the cfg to see if any of them have a sub-cfg that
% would be appropriate to include as a branch
if isempty(datacfg)
  fn = {};
else
  fn = fieldnames(datacfg);
end
for i=1:length(fn)
  if isa(datacfg.(fn{i}), 'struct') && isfield(datacfg.(fn{i}), 'cfg')
    % increment the depth counter
    depth = depth + 1;
    
    % increment the branch counter
    branch = branch + extraBranchesFound;
    extraBranchesFound = 1;
    
    ft_analysisprotocol(cfg, datacfg.(fn{i}).cfg);
  end
end

% revert to the orignal parent
parent = prev;

if depth==1
  % the recursion has finished, we are again at the top level
  
  % record total processing time and maximum memory requirement
  totalproctime = 0;
  maxmemreq = 0;
  
  % the parents were determined while climbing up the tree
  % now it is time to descend and determine the children
  for branch=1:size(info,1)
    for depth=1:size(info,2)
      if ~isempty(info(branch, depth).parent)
        parentbranch = info(branch, depth).parent(1);
        parentdepth  = info(branch, depth).parent(2);
        info(parentbranch, parentdepth).children{end+1} = [branch depth];
        
        if isfield(info(branch,depth).cfg, 'callinfo') && isfield(info(branch,depth).cfg.callinfo, 'proctime')
          totalproctime = totalproctime + info(branch,depth).cfg.callinfo.proctime;
        end
        if isfield(info(branch,depth).cfg, 'callinfo') && isfield(info(branch,depth).cfg.callinfo, 'procmem')
          maxmemreq = max(maxmemreq, info(branch,depth).cfg.callinfo.procmem);
        end
        
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
      
      
      % this seems pointless?
      %disp(commandline);
      
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
    set(gca,'Units','normalized'); % use normalized units
    for branch=1:size(info,1)
      for depth=1:size(info,2)
        plotinfo(cfg,info(branch,depth),size(info,1),size(info,2));
      end
    end
    axis off;
    axis tight;
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
  
  % give a report on the total time used and max. memory required
  if totalproctime > 3600
    proclabel = sprintf('%5.2g hours', totalproctime./3600);
  else
    proclabel = sprintf('%d seconds', round(totalproctime));
  end
  fprintf('\nthe entire analysis pipeline took %s to run\n', proclabel);
  fprintf('the maximum memory requirement of the analysis pipeline was %d MB\n\n',...
    round(maxmemreq./1024./1024));
  
  % clear all persistent variables
  depth  = [];
  branch = [];
  info   = [];
  parent = [];
else
  % this level of recursion has finished, decrease the depth
  depth = depth - 1;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotinfo(cfg,element,numbranch,numdepth)
if isempty(element.name)
  return
end

%  cfg.showinfo   = string or cell array of strings, information to display
%                   in the gui boxes, can be any combination of
%                   'functionname', 'revision', 'matlabversion',
%                   'computername', 'username', 'calltime', 'timeused',
%                   'memused', 'workingdir', 'scriptdir' (default =
%                   'functionname', only display function name)

% create the text information to display
label = {};
for k = 1:numel(cfg.showinfo)
  switch cfg.showinfo{k}
    
    case 'functionname'
      label{end+1} = ['{\bf ' element.name '}'];
      if k == 1 % add blank line if function name is on top, looks nice
        label{end+1} = '';
        firstLabelIsName = 1;
      end
      
    case 'revision'
      if isfield(element.cfg, 'version') && isfield(element.cfg.version, 'id')
        label{end+1} = element.cfg.version.id;
      else
        label{end+1} = '<revision unknown>';
      end
      
    case 'matlabversion'
      if isfield(element.cfg, 'callinfo') && isfield(element.cfg.callinfo, 'matlab')
        label{end+1} = ['MATLAB ' element.cfg.callinfo.matlab];
      else
        label{end+1} = '<MATLAB version unknown>';
      end
      
    case 'computername'
      if isfield(element.cfg, 'callinfo') && isfield(element.cfg.callinfo, 'hostname')
        label{end+1} = ['Computer name: ' element.cfg.callinfo.hostname];
      else
        label{end+1} = '<hostname unknown>';
      end
      
    case 'username'
      if isfield(element.cfg, 'callinfo') && isfield(element.cfg.callinfo, 'user')
        label{end+1} = ['Username: ' element.cfg.callinfo.user];
      else
        label{end+1} = '<username unknown>';
      end
      
    case 'calltime'
      if isfield(element.cfg, 'callinfo') && isfield(element.cfg.callinfo, 'calltime')
        label{end+1} = ['Function called at ' datestr(element.cfg.callinfo.calltime)];
      else
        label{end+1} = '<function call time unknown>';
      end
      
    case 'timeused'
      if isfield(element.cfg, 'callinfo') && isfield(element.cfg.callinfo, 'proctime')
        label{end+1} = sprintf('Function call required %d seconds.', round(element.cfg.callinfo.proctime));
      else
        label{end+1} = '<processing time unknown>';
      end
      
    case 'memused'
      if isfield(element.cfg, 'callinfo') && isfield(element.cfg.callinfo, 'procmem')
        label{end+1} = sprintf('Function call required %d MB.', round(element.cfg.callinfo.procmem/1024/1024));
      else
        label{end+1} = '<memory requirement unknown>';
      end
      
    case 'workingdir'
      if isfield(element.cfg, 'callinfo') && isfield(element.cfg.callinfo, 'pwd')
        label{end+1} = sprintf('Working directory was %s.', element.cfg.callinfo.pwd);
      else
        label{end+1} = '<working directory unknown>';
      end
      
    case 'scriptpath'
      if isfield(element.cfg, 'version') && isfield(element.cfg.version, 'name')
        label{end+1} = sprintf('Full path to script was %s.', element.cfg.version.name);
      else
        label{end+1} = '<script path unknown>';
      end
  end
end

% escape underscores
label = strrep(label, '_', '\_');

% strip blank line if present and not needed
if strcmp(label{end},'')
  label(end) = [];
end

% compute width and height of each box
% note that axis Units are set to Normalized
wh = [1./numbranch 1./numdepth];
boxpadding = wh ./ 5 ./ numel(label);
boxmargin = wh ./ 5;
% adjust actual, inner, width/height for computed padding and margin
wh = wh - boxpadding.*2 - boxmargin.*2;

% create the 4 corners for our patch
x = [0 1 1 0] .* (wh(1)+boxpadding(1).*2);
y = [0 0 1 1] .* (wh(2)+boxpadding(2).*2);
% close the patch
x = [x x(end)];
y = [y y(end)];
% move the patch
location = (element.this-1).*(wh+boxpadding.*2) + element.this.*boxmargin;
% location is at bottom left corner of patch
x = x + location(1);
y = y + location(2);

p = patch(x', y', 0);
set(p, 'Facecolor', [1 1 0.6])

% store data for this patch
tmpGuidata = guidata(p);
if ~isfield(tmpGuidata, 'patches')
  tmpGuidata.patches = {};
end
tmpGuidata.patches{end+1} = [];
tmpGuidata.patches{end}.x = x';
tmpGuidata.patches{end}.y = y';
tmpGuidata.patches{end}.element = element;
guidata(p, tmpGuidata);

if numel(label) == 1
  % center of patch
  textloc = location+boxpadding+wh./2;
  l = text(textloc(1), textloc(2), label);
  set(l, 'HorizontalAlignment', 'center');
  set(l, 'VerticalAlignment', 'middle');
  set(l, 'fontUnits', 'points');
  set(l, 'fontSize', 10);
else
  % top left corner of patch, inside padding
  textloc = [location(1)+boxpadding(1) location(2)+wh(2)];
  l = text(textloc(1), textloc(2), label);
  set(l, 'HorizontalAlignment', 'left');
  set(l, 'VerticalAlignment', 'top');
  set(l, 'fontUnits', 'points');
  set(l, 'fontSize', 10);
end

set(l, 'interpreter', 'tex');

% draw an arrow if appropriate
if ~isempty(element.parent)
  parentlocation = (element.parent-1).*(wh+boxpadding.*2) + element.parent.*boxmargin;
  tip = parentlocation + [0.5 1].*wh + [1 2].*boxpadding;
  base = location + [0.5 0].*wh + [1 0].*boxpadding;
  arrow(base,tip,'Length',8);
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
patches = guidata(h);
patches = patches.patches; % stupid matlab syntax doesn't allow guidata(h).patches

for k = 1:numel(patches)
  patchX = patches{k}.x;
  patchY = patches{k}.y;
  
  if (x >= patchX(1) && x <= patchX(2) ...
      && y >= patchY(1) && y <= patchY(3))
    uidisplaytext(patches{k}.element.script, patches{k}.element.name);
    break;
  end
end

