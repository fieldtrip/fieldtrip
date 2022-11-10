function [neighbours, cfg] = ft_prepare_neighbours(cfg, data)

% FT_PREPARE_NEIGHBOURS finds the channel neighbours for spatial clustering
% or interpolation of bad channels. Using the 'distance' method, neighbours
% are based on a minimum neighbourhood distance (in cfg.neighbourdist).
% Using the 'triangulation' method calculates a triangulation based on a 2D
% projection of the sensor positions. The 'template' method loads a default
% template for the given data type. Alternatively, using the 'parcellation'
% method, in combination with an atlas as input data, spatial neighbours
% of parcels are determined, based on the spatial relationship between the
% labeled mesh vertices. Currently, only atlases defined on a triangular
% mesh are supported.
%
% Use as
%   neighbours = ft_prepare_neighbours(cfg)
% or
%   neighbours = ft_prepare_neighbours(cfg, data)
% with an input data structure with the channels of interest and that
% contains a sensor description, or represents an atlas, see FT_READ_ATLAS
%
% The configuration can contain
%   cfg.channel       = channels in the data for which neighbours should be determined
%   cfg.method        = 'distance', 'triangulation' or 'template'
%   cfg.template      = name of the template file, e.g. CTF275_neighb.mat
%   cfg.neighbourdist = number, maximum distance between neighbouring sensors
%                       (only for 'distance', default is 40 mm)
%   cfg.compress      = 'yes' or 'no', add extra edges by compressing in the
%                       x- and y-direction (only for 'triangulation', default is yes)
%   cfg.feedback      = 'yes' or 'no' (default = 'no')
%
% The 3D sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad          = structure with gradiometer definition or filename, see FT_READ_SENS
%
% The 2D channel positions can be specified as
%   cfg.layout        = filename of the layout, see FT_PREPARE_LAYOUT
%
% With an atlas in the input, the method 'parcellation' has the additional
% options
%   cfg.parcellation  = string that denotes the field in the atlas that is to be used
%
% The output is an array of structures with the "neighbours" which is
% structured like this:
%        neighbours(1).label = 'Fz';
%        neighbours(1).neighblabel = {'Cz', 'F3', 'F3A', 'FzA', 'F4A', 'F4'};
%        neighbours(2).label = 'Cz';
%        neighbours(2).neighblabel = {'Fz', 'F4', 'RT', 'RTP', 'P4', 'Pz', 'P3', 'LTP', 'LT', 'F3'};
%        neighbours(3).label = 'Pz';
%        neighbours(3).neighblabel = {'Cz', 'P4', 'P4P', 'Oz', 'P3P', 'P3'};
%        etc.
%
% Note that a channel is not considered to be a neighbour of itself.
%
% See also FT_NEIGHBOURPLOT, FT_PREPARE_LAYOUT, FT_DATATYPE_SENS,
% FT_READ_SENS, FT_READ_ATLAS

% Copyright (C) 2006-2021, Eric Maris, Jorn M. Horschig, Robert Oostenveld,
% Jan-Mathijs Schoffelen
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

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    data
ft_preamble provenance data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% the data can be passed as input arguments or can be read from disk
hasdata  = exist('data', 'var');
hasatlas = hasdata && (ft_datatype(data, 'mesh+label') || ft_datatype(data, 'volume+label') || ft_datatype(data, 'source+label'));

% these undocumented methods are needed to support some of the high-level FT functions that call this function
if ~isfield(cfg, 'method')
  if ~isfield(cfg, 'neighbours')
    ft_error('cannot figure out how to construct neighbours, please specify cfg.neighbours or cfg.method and call this function directly')
  else
    if ischar(cfg.neighbours)
      ft_notice('reading neighbours from file %s', cfg.neighbours);
      cfg.neighbours = loadvar(cfg.neighbours);
    elseif isstruct(cfg.neighbours) && ~isempty(cfg.neighbours)
      ft_notice('using specified neighbours for the channels');
    elseif isempty(cfg.neighbours) && ~isfield(cfg, 'neighbourdist')
      ft_notice('not using neighbours for the channels');
      % make an empty neighbour structure
      tmp = struct('label', [], 'neighblabel', []);
      cfg.neighbours = tmp([]);
    end
    % this is a hack to get past the case-statement
    cfg.method = 'specified';
  end
end

if hasatlas && ~strcmp(cfg.method, 'parcellation')
  ft_error('an atlas as input data argument requires cfg.method = ''parcellation''');
elseif strcmp(cfg.method, 'parcellation') && ~hasatlas
  ft_error('cfg.method = ''parcellation'' requires an atlas as input data');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'required',   {'method'});
cfg = ft_checkconfig(cfg, 'renamed',    {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed',    {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed',    {'optofile', 'opto'});
cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'tri', 'triangulation'});

% set the defaults
cfg.feedback     = ft_getopt(cfg, 'feedback',     'no');
cfg.channel      = ft_getopt(cfg, 'channel',      'all');
cfg.compress     = ft_getopt(cfg, 'compress',     'yes');
cfg.parcellation = ft_getopt(cfg, 'parcellation', 'parcellation');

if hasdata && ~hasatlas
  % check if the input data is valid for this function
  data = ft_checkdata(data);
  % set the default for senstype depending on the data
  if isfield(data, 'grad')
    cfg.senstype = ft_getopt(cfg, 'senstype', 'meg');
  elseif isfield(data, 'elec')
    cfg.senstype = ft_getopt(cfg, 'senstype', 'eeg');
  elseif isfield(data, 'opto')
    cfg.senstype = ft_getopt(cfg, 'senstype', 'opto');
  else
    cfg.senstype = ft_getopt(cfg, 'senstype', []);
  end
end

if strcmp(cfg.method, 'distance') || strcmp(cfg.method, 'triangulation')
  % these methods require channel positions in either 3D or in 2D
  
  if isfield(cfg, 'layout') && ~isempty(cfg.layout)
    % get 2D positions from the layout
    tmpcfg  = keepfields(cfg, {'layout', 'channel',   'rows', 'columns', 'commentpos', ...
               'skipcomnt', 'scalepos',  'skipscale', 'projection', 'viewpoint', ...
               'rotate',    'width',     'height',    'elec', 'grad', 'opto', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
    tmpcfg.skipscale = 'yes';
    tmpcfg.skipcomnt = 'yes';
    layout  = ft_prepare_layout(tmpcfg);
    chanpos = layout.pos;
    label   = layout.label;
  else
    % get 3D positions from the sensor description
    if hasdata
      sens = ft_fetch_sens(cfg, data);
    else
      sens = ft_fetch_sens(cfg);
    end
    chanpos = sens.chanpos;
    label   = sens.label;
  end
  
  if hasdata
    % remove channels that are not in data
    [dum, sensidx] = match_str(data.label, label);
    chanpos = chanpos(sensidx, :);
    label   = label(sensidx);
  end
  
  % select the desired channels
  desired = ft_channelselection(cfg.channel, label);
  [sensidx] = match_str(label, desired);
  chanpos = chanpos(sensidx, :);
  label   = label(sensidx);
  
  if strcmp(ft_senstype(label), 'neuromag306')
    ft_warning('Please be aware of the different sensor types in neuromag306 system, see http://www.fieldtriptoolbox.org/faq/why_are_there_multiple_neighbour_templates_for_the_neuromag306_system');
  end
  
end % if distance or triangulation

switch cfg.method
  case 'specified'
    % use the neighbours as specified by the user
    neighbours = cfg.neighbours;
    
  case 'template'
    fprintf('Trying to load sensor neighbours from a template\n');
    % determine from where to load the neighbour template
    if ~isfield(cfg, 'template')
      % if data has been put in, try to estimate the sensor type
      if hasdata
        fprintf('Estimating sensor type of data to determine the layout filename\n');
        senstype = ft_senstype(data.label);
        fprintf('Data is of sensor type ''%s''\n', senstype);
        if ~exist([senstype '_neighb.mat'], 'file')
          if exist([senstype '.lay'], 'file')
            cfg.layout = [senstype '.lay'];
          else
            fprintf('Name of sensor type does not match name of layout- and template-file\n');
          end
        else
          cfg.template = [senstype '_neighb.mat'];
        end
      end
    end
    % if that failed
    if ~isfield(cfg, 'template')
      % check whether a layout can be used
      if ~isfield(cfg, 'layout')
        % error if that fails as well
        ft_error('You need to define a template or layout or give data as an input argument when ft_prepare_neighbours is called with cfg.method=''template''');
      end
      fprintf('Using the 2-D layout filename to determine the template filename\n');
      cfg.template = [strtok(cfg.layout, '.') '_neighb.mat'];
    end
    % adjust filename
    if ~exist(cfg.template, 'file')
      cfg.template = lower(cfg.template);
    end
    % add necessary extensions
    if numel(cfg.template) < 4 || ~isequal(cfg.template(end-3:end), '.mat')
      if numel(cfg.template) < 7 || ~isequal(cfg.template(end-6:end), '_neighb')
        cfg.template = [cfg.template, '_neighb'];
      end
      cfg.template = [cfg.template, '.mat'];
    end
    % check for existence
    if ~exist(cfg.template, 'file')
      ft_error('Template file could not be found - please check spelling or see http://www.fieldtriptoolbox.org/faq/how_can_i_define_my_own_neighbourhood_template (please consider sharing it with others via the FT mailing list)');
    end
    load(cfg.template);
    fprintf('Successfully loaded neighbour structure from %s\n', cfg.template);
    
  case 'distance'
    % use a smart default for the distance
    if ~isfield(cfg, 'neighbourdist')
      if exist('sens', 'var')
        sens = ft_determine_units(sens);
        cfg.neighbourdist = 40 * ft_scalingfactor('mm', sens.unit);
      elseif exist('layout', 'var')
        layout = ft_determine_units(layout);
        cfg.neighbourdist = 40 * ft_scalingfactor('mm', layout.unit);
      end
      fprintf('using a distance threshold of %g\n', cfg.neighbourdist);
    end
    neighbours = compneighbstructfrompos(chanpos, label, cfg.neighbourdist);
    
  case 'triangulation'
    if size(chanpos, 2)==2 || all(chanpos(:,3)==0)
      % the sensor positions are already on a 2D plane
      prj = chanpos(:,1:2);
    else
      % project sensor on a 2D plane
      prj = elproj(chanpos);
    end
    % make a 2d delaunay triangulation of the projected points
    tri = delaunay(prj(:,1), prj(:,2));
    if strcmp(cfg.compress, 'yes')
      tri_x = delaunay(prj(:,1)./2, prj(:,2)); % compress in the x-direction
      tri_y = delaunay(prj(:,1), prj(:,2)./2); % compress in the y-direction
      tri = [tri; tri_x; tri_y];
    end
    neighbours = compneighbstructfromtri(chanpos, label, tri);
    
  case 'parcellation'
    
    if ~isfield(data, cfg.parcellation)
      ft_error('required field %s not present in the input atlas', cfg.parcellation);
    end
    
    switch ft_datatype(data)
      case 'mesh+label'
        
        tri = data.tri;
        
        % ensure that the vertices are indexed starting from 1
        if min(tri(:))==0
          tri = tri + 1;
        end
        
        % ensure that the vertices are indexed according to 1:number of unique vertices
        tri = tri_reindex(tri);
        
        % create the unique edges from the triangulation
        edges  = [tri(:,[1 2]); tri(:,[1 3]); tri(:,[2 3])];
        edges  = double(unique(sort([edges; edges(:,[2 1])],2), 'rows'));
        
        % get the parcel values for the edges that 'go across parcels'
        boundary = data.(cfg.parcellation)(edges);
        boundary = boundary(boundary(:,1)~=boundary(:,2),:);
        
        % unique parcel-crossing edges
        uboundary = unique(boundary, 'rows');
        %uboundary = unique([boundary; boundary(:,[2 1])], 'rows');
        
        % count
        nboundary = zeros(size(uboundary, 1), 1);
        for k = 1:numel(nboundary)
          nboundary(k) = sum(boundary(:,1)==uboundary(k,1) & boundary(:,2)==uboundary(k,2));
        end
        
        nthresh   = 2; % at least nthresh edges need to be there, in order to count as a neighbour, e.g.: touching at the corner does not count
        uboundary = uboundary(nboundary>nthresh, :);
        uboundary = [uboundary; uboundary(:,[2 1])]; % to make the adjacency matrix symmetric
        
        % fill the adjacency matrix
        n   = size(uboundary,1);
        adj = full(sparse(uboundary(:,1),uboundary(:,2),ones(n,1)));
        
        label = data.([cfg.parcellation 'label']);
        
        neighbours = struct([]);
        for k = 1:size(adj,1)
          neighbours(k).label       = label{k};
          neighbours(k).neighblabel = label(adj(k,:)>0);
        end
        
      case {'volume+label' 'source+label'}
        ft_error('not yet implemented');
      otherwise
        % this shouldn't happen, but throw an error just to be sure
        ft_error('the input data should be an atlas with cfg.method = ''parcellation''');
    end
    
  otherwise
    ft_error('unsupported method ''%s''', cfg.method);
end

% only select those channels that are in the data
if isfield(cfg, 'channel') && ~isempty(cfg.channel)
  if hasdata && ~hasatlas
    desired = ft_channelselection(cfg.channel, data.label);
  else
    desired = ft_channelselection(cfg.channel, {neighbours(:).label});
  end
elseif hasdata
  desired = data.label;
else
  desired = {neighbours(:).label};
end

if ~isempty(desired)
  complete = struct;
  for i = 1:numel(desired)
    complete(i).label = desired{i};
    sel = find(strcmp({neighbours(:).label}, desired{i}));
    if numel(sel)==1
      % take the set of neighbours from the definition
      complete(i).neighblabel = neighbours(sel).neighblabel;
    else
      % there are no neighbours defined for this channel
      complete(i).neighblabel = {};
    end
  end
  neighbours = complete;
end

for i = 1:length(neighbours)
  % convert them into row-arrays for a nicer code representation with PRINTRSTRUCT
  neighbours(i).neighblabel = neighbours(i).neighblabel(:)';
end
  
k = 0;
for i = 1:length(neighbours)
  if isempty(neighbours(i).neighblabel)
    ft_warning('no neighbours found for %s', neighbours(i).label);
  end
  k = k + length(neighbours(i).neighblabel);
end

if k==0
  ft_warning('No neighbouring channels were specified or found');
else
  fprintf('there are on average %.1f neighbours per channel\n', k/length(neighbours));
end

if strcmp(cfg.feedback, 'yes') && ~hasatlas
  % give some graphical feedback
  tmpcfg = keepfields(cfg, {'layout', 'rows', 'columns', 'commentpos', 'skipcomnt', 'scalepos', 'skipscale', 'projection', 'viewpoint', 'rotate', 'width', 'height', 'elec', 'grad', 'opto', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
  tmpcfg.neighbours = neighbours;
  if hasdata
    tmpcfg.senstype = cfg.senstype;
    ft_neighbourplot(tmpcfg, data);
  else
    ft_neighbourplot(tmpcfg);
  end
elseif strcmp(cfg.feedback, 'yes')
  ft_notice('no visual feedback about the neighbours is possible with the current input');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   data
ft_postamble provenance neighbours
ft_postamble history    neighbours
ft_postamble savevar    neighbours


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that compute the neighbourhood geometry from the
% gradiometer/electrode positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neighbours] = compneighbstructfrompos(chanpos, label, neighbourdist)

nchan = length(label);

% compute the distance between all sensors
dist = zeros(nchan,nchan);
for i = 1:nchan
  dist(i,:) = sqrt(sum((chanpos(1:nchan,:) - repmat(chanpos(i,:), nchan, 1)).^2,2))';
end

% find the neighbouring electrodes based on distance
% later we have to restrict the neighbouring electrodes to those actually selected in the dataset
channeighbstructmat = (dist<neighbourdist);

% electrode istelf is not a neighbour
channeighbstructmat = (channeighbstructmat .* ~eye(nchan));

% convert back to logical
channeighbstructmat = logical(channeighbstructmat);

% construct a struct-array with all neighbours
neighbours=struct;
for i = 1:nchan
  neighbours(i).label       = label{i};
  neighbours(i).neighblabel = label(channeighbstructmat(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that computes the neighbourhood geometry from the
% triangulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neighbours] = compneighbstructfromtri(chanpos, label, tri)

nchan = length(label);

channeighbstructmat = zeros(nchan,nchan);
% mark neighbours according to triangulation
for i=1:size(tri, 1)
  channeighbstructmat(tri(i, 1), tri(i, 2)) = 1;
  channeighbstructmat(tri(i, 1), tri(i, 3)) = 1;
  channeighbstructmat(tri(i, 2), tri(i, 1)) = 1;
  channeighbstructmat(tri(i, 3), tri(i, 1)) = 1;
  channeighbstructmat(tri(i, 2), tri(i, 3)) = 1;
  channeighbstructmat(tri(i, 3), tri(i, 2)) = 1;
end

% construct a structured cell-array with all neighbours
neighbours = struct;
alldist = [];
for i = 1:nchan
  neighbours(i).label       = label{i};
  neighbidx                 = find(channeighbstructmat(i,:));
  neighbours(i).dist        = sqrt(sum((repmat(chanpos(i, :), numel(neighbidx), 1) - chanpos(neighbidx, :)).^2, 2));
  alldist                   = [alldist; neighbours(i).dist];
  neighbours(i).neighblabel = label(neighbidx);
end

% remove neighbouring channels that are too far away (imporntant e.g. in
% case of missing sensors)
neighbdist = mean(alldist)+3*std(alldist);
for i=1:nchan
  idx = neighbours(i).dist > neighbdist;
  neighbours(i).dist(idx)         = [];
  neighbours(i).neighblabel(idx)  = [];
end
neighbours = rmfield(neighbours, 'dist');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that ensures a consistent triangulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newtri] = tri_reindex(tri)

% this subfunction reindexes tri such that they run from 1:number of unique vertices
newtri       = tri;
[srt, indx]  = sort(tri(:));
tmp          = cumsum(double(diff([0;srt])>0));
newtri(indx) = tmp;
