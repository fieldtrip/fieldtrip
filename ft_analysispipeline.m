function [pipeline] = ft_analysispipeline(cfg, data)

% FT_ANALYSIPIPELINE reconstructs the complete analysis pipeline that was used to create
% the input FieldTrip data structure. The pipeline will be visualized as a flowchart.
% In the future it will be possible to output the complete pipeline as a MATLAB script
% or in a specialized pipeline format (e.g. PSOM, JIST, LONI, Taverna).
%
% Use as
%   output = ft_analysispipeline(cfg, data)
%
% The first cfg input contains the settings that apply to the behaviour of this
% particular function and the second data input argument can be the output of any
% FieldTrip function, e.g. FT_PREPROCESSING, FT_TIMELOCKANALYSIS, FT_SOURCEANALYSIS,
% FT_FREQSTATISTICS or whatever you like.
%
% Alternatively, for the second input argument you can also only give the configuration
% of the processed data (i.e. "data.cfg") instead of the full data.
%
% The configuration options that apply to the behaviour of this function are
%   cfg.filename   = string, filename without the extension
%   cfg.filetype   = string, can be 'matlab', 'html' or 'dot'
%   cfg.feedback   = string, 'no', 'text', 'gui' or 'yes', whether text and/or
%                    graphical feedback should be presented (default = 'yes')
%   cfg.showinfo   = string or cell array of strings, information to display
%                    in the gui boxes, can be any combination of
%                    'functionname', 'revision', 'matlabversion',
%                    'computername', 'username', 'calltime', 'timeused',
%                    'memused', 'workingdir', 'scriptpath' (default =
%                    'functionname', only display function name). Can also
%                    be 'all', show all pipeline. Please note that if you want
%                    to show a lot of information, this will require a lot
%                    of screen real estate.
%   cfg.remove     = cell-array with strings, determines which objects will
%                    be removed from the configuration prior to writing it to
%                    file. For readibility of the script, you may want to
%                    remove the large objectssuch as event structure, trial
%                    definition, source positions
%  cfg.keepremoved = 'yes' or 'no', determines whether removed fields are
%                    completely removed, or only replaced by a short textual
%                    description (default = 'no')
%
% This function uses the nested cfg and cfg.previous that are present in
% the data structure. It will use the configuration and the nested previous
% configurations to climb all the way back into the tree. This funtction
% will print a complete MATLAB script to screen (and optionally to file).
% Furthermore, it will show an interactive graphical flowchart
% representation of the steps taken during the pipeline(i). In the flowchart
% you can click on one of the steps to see the configuration details of
% that pipeline(i).
%
% Note that the nested cfg and cfg.previous in your data might not contain
% all details that are required to reconstruct a complete and valid
% analysis script.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this, the input data will be read from a *.mat file on disk. The
% file should contain only a single variable, corresponding with the input structure.
%
% See also FT_PREPROCESSING, FT_TIMELOCKANALYSIS, FT_FREQANALYSIS, FT_SOURCEANALYSIS,
% FT_CONNECTIVITYANALYSIS, FT_NETWORKANALYSIS

% Copyright (C) 2014-2015, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% callinfo feedback is highly annoying in this recursive function
% do this here, otherwise ft_defaults will override our setting
if ~isfield(cfg, 'showcallinfo'), cfg.showcallinfo = 'no';   end

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the defaults
cfg.filename    = ft_getopt(cfg, 'filename');
cfg.showinfo    = ft_getopt(cfg, 'showinfo', {'functionname'});
cfg.keepremoved = ft_getopt(cfg, 'keepremoved', 'no');
cfg.feedback    = ft_getopt(cfg, 'feedback', 'text');
cfg.prune       = ft_getopt(cfg, 'prune', 'yes');
cfg.filetype    = ft_getopt(cfg, 'filetype');
cfg.fontsize    = ft_getopt(cfg, 'fontsize', 10);


if isempty(cfg.filetype) && ~isempty(cfg.filename)
  [p, f, x] = fileparts(cfg.filename);
  switch x
    case '.m'
      cfg.filetype = 'matlab';
    case '.html'
      cfg.filetype = 'html';
    case '.dot'
      cfg.filetype = 'dot';
    otherwise
      error('cannot determine filetype');
  end
end

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
    'vol.bnd.pos'
    'vol.bnd.tri'
    'headmodel.bnd.pos'
    'headmodel.bnd.tri'
    };
elseif ~iscell(cfg.remove)
  cfg.remove = {cfg.remove};
end

if strcmp(cfg.showinfo, 'all')
  cfg.showinfo = {
    'functionname'
    'revision'
    'matlabversion'
    'computername'
    'architecture'
    'username'
    'calltime'
    'timeused'
    'memused'
    'workingdir'
    'scriptpath'
    };
end

if ~isfield(cfg, 'showinfo')
  cfg.showinfo = {'functionname'};
elseif ~iscell(cfg.showinfo)
  cfg.showinfo = {cfg.showinfo};
end

% we are only interested in the cfg-part of the data
if isfield(data, 'cfg')
  datacfg = data.cfg;
else
  datacfg = data;
end
clear data

% walk the tree, gather information about each node
ft_progress('init', cfg.feedback, 'parsing provenance...');
pipeline = walktree(datacfg);
ft_progress('close');

% convert the cell array into a structure array
for i=1:length(pipeline)
  tmp(i) = pipeline{i};
end
pipeline = tmp;

if istrue(cfg.prune)
  % prune the double occurences
  [dummy, indx] = unique({pipeline.this});
  pipeline = pipeline(sort(indx));
end

% start at the end of the tree and determine the level of each of the parents
hasparent = false(size(pipeline));
haschild  = false(size(pipeline));
for i=1:length(pipeline)
  hasparent(i) = ~isempty(pipeline(i).parent);
  haschild(i)  = any(strcmp(pipeline(i).this, [pipeline.parent]));
end

% construct a matrix with all pipeline steps
width  = zeros(size(pipeline));
height = zeros(size(pipeline));
level  = 1;
% the items without children start at height 1
sel = find(~haschild);
while ~isempty(sel)
  height(sel) = level;
  % find the parents of the items at this level
  sel = match_str({pipeline.this}, [pipeline(sel).parent]);
  % the parents should be at least one level higher
  height(sel) = level + 1;
  % continue with the next level
  level = level + 1;
end
for i=1:max(height)
  sel = find(height==i);
  width(sel) = 1:length(sel);
end
for i=1:length(pipeline)
  pipeline(i).position = [height(i) width(i)];
end

% sort according to a decreasing level, i.e. the last pipeline(i) at the end
[dummy, indx] = sortrows(-[height(:) width(:)]);
pipeline = pipeline(indx);

if isempty(cfg.filename)
  pipeline2matlabfigure(cfg, pipeline);
else
  switch cfg.filetype
    case 'matlab'
      pipeline2matlabscript(cfg, pipeline);
    case 'dot'
      pipeline2dotfile(cfg, pipeline);
    case 'html'
      pipeline2htmlfile(cfg, pipeline);
    otherwise
      error('unsupported filetype');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for recursive walking along the cfg.previous.previous info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function info = walktree(cfg)

if isempty(cfg) || (isstruct(cfg) && numel(fieldnames(cfg))==0)
  info = [];
  return
end

this = getnode(cfg);

% parse all previous steps
if isfield(cfg, 'previous') && ~isempty(cfg.previous) && iscell(cfg.previous)
  previous = cellfun(@walktree, cfg.previous, 'UniformOutput', false);
  if iscell(previous{1})
    previous = cat(2, previous{:});
  end
elseif isfield(cfg, 'previous') && ~isempty(cfg.previous) && isstruct(cfg.previous)
  previous = walktree(cfg.previous);
elseif isfield(cfg, 'previous') && ~isempty(cfg.previous)
  error('unexpected content in cfg.previous');
else
  previous = {};
end

% parse the side branches, e.g. cfg.vol/cfg.headmodel and cfg.layout
fn = fieldnames(cfg);
branch = {};
for i=1:numel(fn)
  if isstruct(cfg.(fn{i})) && isfield(cfg.(fn{i}), 'cfg')
    branch = [walktree(cfg.(fn{i}).cfg) branch];
    this.parent{end+1} = branch{1}.this; % the start of the branch is a parent to this element
  end
end

ft_progress(rand(1), 'parsing provenance for %s\n', this.name); % FIXME no percentage complete known
drawnow

% the order of the output elements matters for the recursion
info = [{this} branch previous];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for gathering the information about each pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function node = getnode(cfg)
[p, f, x]      = myfileparts(getvalue(cfg, 'version.name'));
node.cfg       = cfg;
node.name      = f;
node.id        = getvalue(cfg, 'version.id');
node.this      = ft_hash(cfg);
if isfield(cfg, 'previous') && ~isempty(cfg.previous) && iscell(cfg.previous)
  % skip the entries that are empty
  cfg.previous = cfg.previous(~cellfun(@isempty, cfg.previous));
  node.parent   = cellfun(@ft_hash, cfg.previous, 'UniformOutput', false);
elseif isfield(cfg, 'previous') && ~isempty(cfg.previous) && isstruct(cfg.previous)
  node.parent   = {ft_hash(cfg.previous)};
elseif isfield(cfg, 'previous') && ~isempty(cfg.previous)
  error('unexpected content in cfg.previous');
else
  node.parent   = {};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = getvalue(s, f)
if issubfield(s, f)
  v = getsubfield(s, f);
else
  v = 'unknown';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, f, x] = myfileparts(filename)
if ispc
  filename(filename=='/') = filesep;
else
  filename(filename=='\') = filesep;
end
[p, f, x] = fileparts(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = button(h, eventdata, handles, varargin)
pos = get(get(gcbo, 'CurrentAxes'), 'CurrentPoint');
x = pos(1,1);
y = pos(1,2);
pipeline = guidata(h);

for i=1:numel(pipeline)
  if (x >= pipeline(i).x(1) && x <= pipeline(i).x(2) && y >= pipeline(i).y(1) && y <= pipeline(i).y(3))
    cfg = pipeline(i).cfg;
    if isfield(cfg, 'previous')
      cfg = rmfield(cfg, 'previous');
    end
    % use a helper function to remove uninteresting fields
    cfg = removefields(cfg, ignorefields('pipeline'), 'recursive', 'yes');
    % use a helper function to remove too large fields
    cfg.checksize = 3000;
    cfg = ft_checkconfig(cfg, 'checksize', 'yes');
    cfg = rmfield(cfg, 'checksize');
    script = printstruct('cfg', cfg);
    uidisplaytext(script, pipeline(i).name);
    break;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pipeline2matlabfigure(cfg, pipeline)

fprintf('plotting pipeline as MATLAB figure\n');

layout = cell(numel(pipeline));
for i=1:length(pipeline)
  % get the vertical and horizontal position in integer values
  % high numbers are towards the begin, low numbers are towards the end of the pipeline
  v = pipeline(i).position(1);
  h = pipeline(i).position(2);
  layout{v,h} = pipeline(i).this;
end

% remove the empty columns
maxheight = max(sum(~cellfun(@isempty, layout),1));
maxwidth  = max(sum(~cellfun(@isempty, layout),2));
layout    = layout(1:maxheight,1:maxwidth);

fig = figure;
hold on
axis manual; % the axis should not change during the contruction of the arrows, otherwise the arrowheads will be distorted
set(gca,'Units','normalized'); % use normalized units
set(gcf, 'ToolBar', 'none');
axis([0 1 0 1])
axis off;
axis tight;

for i=1:numel(pipeline)
  
  label = makelabel(pipeline(i), cfg.showinfo);
  
  % dublicate backslashes to escape tex interpreter (in case of windows filenames)
  label = strrep(label, '\', '\\');
  label = strrep(label, '{\\bf', '{\bf'); % undo for bold formatting
  
  % escape underscores
  label = strrep(label, '_', '\_');
  
  % strip blank line if present and not needed
  if strcmp(label{end},'')
    label(end) = [];
  end
  
  % compute width and height of each box, note that axis Units are set to Normalized
  boxsize = 1./[maxwidth+1 maxheight+3];
  
  % create the 4 corners for our patch, close the patch by returning to the start point
  x = ([0 1 1 0 0]-0.5) .* boxsize(1);
  y = ([0 0 1 1 0]-0.5) .* boxsize(2);
  
  % position the patch
  location    = pipeline(i).position([2 1]);
  location(1) = (location(1)-0.5)/maxwidth;
  location(2) = (location(2)-0.5)/maxheight;
  
  % the location specifies the center of the patch
  x = x + location(1);
  y = y + location(2);
  
  p = patch(x', y', 0);
  set(p, 'Facecolor', [1 1 0.6])
  
  pipeline(i).x = x;
  pipeline(i).y = y;
  guidata(fig, pipeline);
  
  if length(label)==1
    textloc = location;
    l = text(textloc(1), textloc(2), label);
    set(l, 'HorizontalAlignment', 'center');
    set(l, 'VerticalAlignment', 'middle');
    set(l, 'fontUnits', 'points');
    set(l, 'fontSize', cfg.fontsize);
    set(l, 'interpreter', 'tex');
  else
    textloc = location;
    textloc(1) = textloc(1)-boxsize(1)/2;
    textloc(2) = textloc(2)+boxsize(2)/2;
    
    l = text(textloc(1), textloc(2), label);
    set(l, 'HorizontalAlignment', 'left');
    set(l, 'VerticalAlignment', 'top');
    set(l, 'fontUnits', 'points');
    set(l, 'fontSize', cfg.fontsize);
    set(l, 'interpreter', 'tex');
  end
  
  % draw an arrow if appropriate
  n = length(pipeline(i).parent);
  
  for j=1:n
    [parentlocation(2), parentlocation(1)] = ind2sub([maxheight, maxwidth], find(strcmp(layout(:), pipeline(i).parent{j}), 1, 'first'));
    % parentlocation = info(find(strcmp({pipeline.this}, analysis.parent{j}), 1, 'first')).position;
    parentlocation(1) = (parentlocation(1)-0.5)/maxwidth;
    parentlocation(2) = (parentlocation(2)-0.5)/maxheight;
    base = parentlocation + [0 -0.5].*boxsize;
    if false % n>1
      % distribute the inputs evenly over the box
      rel = 2*(j+n-1)/(n-1)-1; % relative location between -1.0 and 1.0
      rel = rel*(n-1)/(n+3);   % compress the relative location
      tip = location + [0 0.5].*boxsize + [rel 0].*boxsize/2;
    else
      % put it in the centre
      tip  = location + [0 0.5].*boxsize;
    end
    arrow(base, tip, 'length', 8, 'lineWidth', 1);
  end
  
  
end % for numel(info)

set(fig, 'WindowButtonUpFcn', @button);
% set(fig, 'KeyPressFcn', @key);

% add a context menu to the figure
% ftmenu = uicontextmenu; set(gcf, 'uicontextmenu', ftmenu)

% add a regular menu item to the figure
ftmenu  = uimenu(fig, 'Label', 'FieldTrip');
% ftmenu1 = uimenu(ftmenu, 'Label', 'Save pipeline');
% ftmenu2 = uimenu(ftmenu, 'Label', 'Share pipeline');
uimenu(ftmenu, 'Label', 'About',  'Separator', 'on', 'Callback', @menu_about);
% uimenu(ftmenu1, 'Label', 'Save as MATLAB script');
% uimenu(ftmenu1, 'Label', 'Save as PSOM pipeline');
% uimenu(ftmenu1, 'Label', 'Save as HTML page');
% uimenu(ftmenu2, 'Label', 'Share within DCCN');
% uimenu(ftmenu2, 'Label', 'Share on PasteBin.com');
% uimenu(ftmenu2, 'Label', 'Share on MyExperiment.org');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label = makelabel(pipeline, showinfo)
% create the text information to display
label = {};
for k = 1:numel(showinfo)
  switch showinfo{k}
    
    case 'functionname'
      % label{end+1} = ['{\bf ' pipeline(i).name '}'];
      label{end+1} = pipeline.name;
      if k == 1 % add blank line if function name is on top, looks nice
        label{end+1} = '';
      end
      
    case 'revision'
      if isfield(pipeline.cfg, 'version') && isfield(pipeline.cfg.version, 'id')
        label{end+1} = pipeline.cfg.version.id;
      else
        label{end+1} = '<revision unknown>';
      end
      
    case 'matlabversion'
      if isfield(pipeline.cfg, 'callinfo') && isfield(pipeline.cfg.callinfo, 'matlab')
        label{end+1} = ['MATLAB ' pipeline.cfg.callinfo.matlab];
      else
        label{end+1} = '<MATLAB version unknown>';
      end
      
    case 'computername'
      if isfield(pipeline.cfg, 'callinfo') && isfield(pipeline.cfg.callinfo, 'hostname')
        label{end+1} = ['Hostname: ' pipeline.cfg.callinfo.hostname];
      else
        label{end+1} = '<hostname unknown>';
      end
      
    case 'architecture'
      if isfield(pipeline.cfg, 'callinfo') && isfield(pipeline.cfg.callinfo, 'hostname')
        label{end+1} = ['Architecture: ' pipeline.cfg.callinfo.computer];
      else
        label{end+1} = '<architecture unknown>';
      end
      
    case 'username'
      if isfield(pipeline.cfg, 'callinfo') && isfield(pipeline.cfg.callinfo, 'user')
        label{end+1} = ['Username: ' pipeline.cfg.callinfo.user];
      else
        label{end+1} = '<username unknown>';
      end
      
    case 'calltime'
      if isfield(pipeline.cfg, 'callinfo') && isfield(pipeline.cfg.callinfo, 'calltime')
        label{end+1} = ['Function called at ' datestr(pipeline.cfg.callinfo.calltime)];
      else
        label{end+1} = '<function call time unknown>';
      end
      
    case 'timeused'
      if isfield(pipeline.cfg, 'callinfo') && isfield(pipeline.cfg.callinfo, 'proctime')
        label{end+1} = sprintf('Function call required %d seconds', round(pipeline.cfg.callinfo.proctime));
      else
        label{end+1} = '<processing time unknown>';
      end
      
    case 'memused'
      if isfield(pipeline.cfg, 'callinfo') && isfield(pipeline.cfg.callinfo, 'procmem')
        label{end+1} = sprintf('Function call required %d MB', round(pipeline.cfg.callinfo.procmem/1024/1024));
      else
        label{end+1} = '<memory requirement unknown>';
      end
      
    case 'workingdir'
      if isfield(pipeline.cfg, 'callinfo') && isfield(pipeline.cfg.callinfo, 'pwd')
        label{end+1} = sprintf('Working directory was %s', pipeline.cfg.callinfo.pwd);
      else
        label{end+1} = '<working directory unknown>';
      end
      
    case 'scriptpath'
      if isfield(pipeline.cfg, 'version') && isfield(pipeline.cfg.version, 'name')
        label{end+1} = sprintf('Full path to script was %s', pipeline.cfg.version.name);
      else
        label{end+1} = '<script path unknown>';
      end
  end
end % for numel(showinfo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pipeline2matlabscript(cfg, pipeline)

[p, f, x] = fileparts(cfg.filename);
filename = fullfile(p, [f '.m']);

fprintf('exporting MATLAB script to file ''%s''\n', filename);

varname = {};
varhash = {};

for i=1:length(pipeline)
  varname{end+1} = sprintf('variable_%d_%d', pipeline(i).position(1), pipeline(i).position(2));
  varhash{end+1} = pipeline(i).this;
end

for i=1:length(pipeline)
  pipeline(i).inputvar  = {};
  pipeline(i).outputvar = varname{strcmp(varhash, pipeline(i).this)};
  for j=1:length(pipeline(i).parent)
    pipeline(i).inputvar{j} = varname{strcmp(varhash, pipeline(i).parent{j})};
  end
end

sep    = sprintf('\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
script = sep;
for i=1:length(pipeline)
  thiscfg = pipeline(i).cfg;
  if isfield(thiscfg, 'previous')
    thiscfg = rmfield(thiscfg, 'previous');
  end
  cfgstr    = printstruct('cfg', thiscfg);
  if ~isempty(pipeline(i).inputvar)
    inputvar = sprintf(', %s', pipeline(i).inputvar{:});
  else
    inputvar = '';
  end
  cmd       = sprintf('\n%s = %s(cfg%s);', pipeline(i).outputvar, pipeline(i).name, inputvar);
  script    = cat(2, script, cfgstr, cmd, sep);
end

% write the complete script to file
fid = fopen(filename, 'wb');
fprintf(fid, '%s', script);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pipeline2dotfile(cfg, pipeline)

[p, f, x] = fileparts(cfg.filename);
filename = fullfile(p, [f '.dot']);

fprintf('exporting DOT file to ''%s''\n', filename);

% write the complete script to file
fid = fopen(filename, 'wb');

fprintf(fid, 'digraph {\n');

varhash = {pipeline.this};

for i=1:length(pipeline)
  for j=1:length(pipeline(i).parent)
    fprintf(fid, '%d -> %d\n', find(strcmp(varhash, pipeline(i).parent{j})), i);
  end
end

for i=1:length(pipeline)
  label = makelabel(pipeline(i), cfg.showinfo);
  if numel(label)>2
    % left justified
    label = sprintf('%s\\l', label{:});
    label = label(1:end-2);
  else
    % centered
    label = sprintf('%s\\n', label{:});
    label = label(1:end-2);
  end
  fprintf(fid, '%d [label="%s",shape=box,fontsize=%d,URL="http://www.fieldtriptoolbox.org/reference/%s"]\n', i, label, cfg.fontsize, pipeline(i).name);
end

fprintf(fid, '}\n');

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pipeline2htmlfile(cfg, pipeline)

[p, f, x] = fileparts(cfg.filename);
filename = fullfile(p, [f '.html']);

% skip the data-like fields and the fields that probably were not added by the user himself
skipfields = {'previous', 'grid', 'headmodel', 'event', 'warning', 'progress', 'trackconfig', 'checkconfig', 'checksize', 'showcallinfo', 'debug', 'outputfilepresent', 'trackcallinfo', 'trackdatainfo', 'trackusage'};

fprintf('exporting HTML file to ''%s''\n', filename);

html = '';
totalproctime = 0;

ft_progress('init', cfg.feedback, 'serialising cfg-structures...');

for k = 1:numel(pipeline)
  ft_progress(k/numel(pipeline), 'serialising cfg-structure %d from %d', k, numel(pipeline));
  
  % strip away the cfg.previous fields, and all data-like fields
  tmpcfg = removefields(pipeline(k).cfg, skipfields);
  
  usercfg = [];
  
  % record the usercfg and proctime if present
  if isfield(tmpcfg, 'callinfo')
    if isfield(tmpcfg.callinfo, 'usercfg')
      usercfg = removefields(tmpcfg.callinfo.usercfg, skipfields);
      
      % avoid processing usercfg twice
      tmpcfg.callinfo = rmfield(tmpcfg.callinfo, 'usercfg');
    end
    
    if isfield(tmpcfg.callinfo, 'proctime')
      totalproctime = totalproctime + tmpcfg.callinfo.proctime;
    end
  end
  
  html = [html sprintf('nodes["%s"] = {"id": "%s", "name": "%s", "cfg": "%s", "usercfg": "%s", "parentIds": [',...
    pipeline(k).this, pipeline(k).this, pipeline(k).name, escapestruct(tmpcfg), escapestruct(usercfg))];
  
  if ~isempty(pipeline(k).parent)
    for j = 1:numel(pipeline(k).parent)
      html = [html '"' pipeline(k).parent{j} '"'];
      if j < numel(pipeline(k).parent)
        html = [html ','];
      end
    end
  end
  
  html = [html sprintf(']};\n')];
  
  if k == numel(pipeline)
    % we are at the single leaf node
    html = [html sprintf('var leafId = "%s";\n', pipeline(k).this)];
  end
  
end
ft_progress('close');

html = [html(1:end-2) sprintf('\n')];

% load the skeleton and put in the html code
thispath = fileparts(mfilename('fullpath'));
htmlfile = fileread([thispath '/private/pipeline-skeleton.html']);

timestamp = [datestr(now(), 'ddd') ' ' datestr(now())];
proctimestr = print_tim(totalproctime);
htmlfile = strrep(htmlfile, '${TIMESTAMP}', timestamp);
htmlfile = strrep(htmlfile, '${PIPELINE}', html);
htmlfile = strrep(htmlfile, '${USER}', getusername());
htmlfile = strrep(htmlfile, '${PROCTIME}', proctimestr);

% write the file
fid = fopen(filename, 'w');
fwrite(fid, htmlfile, 'uchar');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cfghtml = escapestruct(tmpcfg)
% convert the cfg structure to a suitable string (escape newlines and
% quotes)

% strip away big numeric fields
if isstruct(tmpcfg)
  fields = fieldnames(tmpcfg);
  for k = 1:numel(fields)
    if isnumeric(tmpcfg.(fields{k})) && numel(tmpcfg.(fields{k})) > 400
      tmpcfg.(fields{k}) = '[numeric data of &gt;400 elements stripped]';
    end
  end
end
cfghtml = strrep(printstruct('cfg', tmpcfg), '\', '\\');
cfghtml = strrep(cfghtml, sprintf('\r'), '\r');
cfghtml = strrep(cfghtml, sprintf('\n'), '\n');
cfghtml = strrep(cfghtml, sprintf('\t'), '\t');
cfghtml = strrep(cfghtml, '"', '\"');
