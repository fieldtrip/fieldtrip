function [neighbours, cfg] = ft_prepare_neighbours(cfg, data)

% FT_PREPARE_NEIGHBOURS finds the neighbours of the channels based on three
% different methods. Using the 'distance'-method, prepare_neighbours is
% based on a minimum neighbourhood distance (in cfg.neighbourdist). The
% 'triangulation'-method calculates a triangulation based on a
% two-dimenstional projection of the sensor position. The 'template'-method
% loads a default template for the given data type. prepare_neighbours
% should be verified using cfg.feedback ='yes' or by calling
% ft_neighbourplot
%
% The positions of the channel are specified in a gradiometer or electrode configuration or
% from a layout. The sensor configuration can be passed into this function in three ways:
%  (1) in a configuration field,
%  (2) in a file whose name is passed in a configuration field, and that can be imported using FT_READ_SENS, or
%  (3) in a data field.
%
% Use as
%   neighbours = ft_prepare_neighbours(cfg, data)
%
% The configuration can contain
%   cfg.method        = 'distance', 'triangulation' or 'template'
%   cfg.neighbourdist = number, maximum distance between neighbouring sensors (only for 'distance')
%   cfg.template      = name of the template file, e.g. CTF275_neighb.mat
%   cfg.layout        = filename of the layout, see FT_PREPARE_LAYOUT
%   cfg.channel       = channels for which neighbours should be found
%   cfg.feedback      = 'yes' or 'no' (default = 'no')
%
% The EEG or MEG sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions, see FT_DATATYPE_SENS
%   cfg.grad          = structure with gradiometer definition, see FT_DATATYPE_SENS
%   cfg.elecfile      = name of file containing the electrode positions, see FT_READ_SENS
%   cfg.gradfile      = name of file containing the gradiometer definition, see FT_READ_SENS
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
% Note that a channel is not considered to be a neighbour of itself.
%
% See also FT_NEIGHBOURPLOT

% Copyright (C) 2006-2011, Eric Maris, Jorn M. Horschig, Robert Oostenveld
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
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', {'method'});

% set the defaults
cfg.feedback = ft_getopt(cfg, 'feedback', 'no');
cfg.channel = ft_getopt(cfg, 'channel', 'all');

% the data can be passed as input arguments or can be read from disk
hasdata = exist('data', 'var');

if hasdata
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

if strcmp(cfg.method, 'template')
  neighbours = [];
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

else
  % get the the grad or elec if not present in the data
  if hasdata
    sens = ft_fetch_sens(cfg, data);
  else
    sens = ft_fetch_sens(cfg);
  end
  
  if strcmp(ft_senstype(sens), 'neuromag306')
    ft_warning('Neuromag306 system detected - be aware of different sensor types, see http://www.fieldtriptoolbox.org/faq/why_are_there_multiple_neighbour_templates_for_the_neuromag306_system');
  end
  chanpos = sens.chanpos;
  label   = sens.label;
  
  if nargin > 1
    % remove channels that are not in data
    [dataidx, sensidx] = match_str(data.label, label);
    chanpos = chanpos(sensidx, :);
    label   = label(sensidx);
  end
  
  if ~strcmp(cfg.channel, 'all')
    desired = ft_channelselection(cfg.channel, label);
    [sensidx] = match_str(label, desired);
    chanpos = chanpos(sensidx, :);
    label   = label(sensidx);
  end
  
  switch lower(cfg.method)
    case 'distance'
      % use a smart default for the distance
      if ~isfield(cfg, 'neighbourdist')
        sens = ft_checkdata(sens, 'hasunit', 'yes');
        cfg.neighbourdist = 40 * ft_scalingfactor('mm', sens.unit);
        fprintf('using a distance threshold of %g\n', cfg.neighbourdist);
      end
      
      neighbours = compneighbstructfromgradelec(chanpos, label, cfg.neighbourdist);
    case {'triangulation', 'tri'} % the latter for reasons of simplicity
      if size(chanpos, 2)==2 || all(chanpos(:,3)==0)
        % the sensor positions are already on a 2D plane
        prj = chanpos(:,1:2);
      else
        % project sensor on a 2D plane
        prj = elproj(chanpos);
      end
      % make a 2d delaunay triangulation of the projected points
      tri = delaunay(prj(:,1), prj(:,2));
      tri_x = delaunay(prj(:,1)./2, prj(:,2));
      tri_y = delaunay(prj(:,1), prj(:,2)./2);
      tri = [tri; tri_x; tri_y];
      neighbours = compneighbstructfromtri(chanpos, label, tri);
    otherwise
      ft_error('Method ''%s'' not known', cfg.method);
  end
end

% removed as from Nov 09 2011 - hope there are no problems with this
% if iscell(neighbours)
%   ft_warning('Neighbourstructure is in old format - converting to structure array');
%   neighbours = fixneighbours(neighbours);
% end

% only select those channels that are in the data
neighb_chans = {neighbours(:).label};
if isfield(cfg, 'channel') && ~isempty(cfg.channel)
  if hasdata
    desired = ft_channelselection(cfg.channel, data.label);
  else
    desired = ft_channelselection(cfg.channel, neighb_chans);
  end
elseif (hasdata)
  desired = data.label;
else
  desired = neighb_chans;
end

% in any case remove SCALE and COMNT
desired = ft_channelselection({'all', '-SCALE', '-COMNT'}, desired);

neighb_idx = ismember(neighb_chans, desired);
neighbours = neighbours(neighb_idx);

k = 0;
for i=1:length(neighbours)
  if isempty(neighbours(i).neighblabel)
    ft_warning('FIELDTRIP:NoNeighboursFound', 'no neighbours found for %s\n', neighbours(i).label);
    % JMH: I removed this in Feb 2013 - this is handled above now
    % note however that in case of using a template, this function behaves
    % differently now (neighbourschans can still be channels not in
    % cfg.channel)
    %else % only selected desired channels
    %  neighbours(i).neighblabel = neighbours(i).neighblabel(ismember(neighbours(i).neighblabel, desired));
  end
  k = k + length(neighbours(i).neighblabel);
end

if k==0
  ft_error('No neighbours were found!');
end

fprintf('there are on average %.1f neighbours per channel\n', k/length(neighbours));

if strcmp(cfg.feedback, 'yes')
  % give some graphical feedback
  cfg.neighbours = neighbours;
  if hasdata
    ft_neighbourplot(cfg, data);
  else
    ft_neighbourplot(cfg);
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance neighbours
ft_postamble history    neighbours


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that compute the neighbourhood geometry from the
% gradiometer/electrode positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neighbours] = compneighbstructfromgradelec(chanpos, label, neighbourdist)

nsensors = length(label);

% compute the distance between all sensors
dist = zeros(nsensors,nsensors);
for i=1:nsensors
  dist(i,:) = sqrt(sum((chanpos(1:nsensors,:) - repmat(chanpos(i,:), nsensors, 1)).^2,2))';
end;

% find the neighbouring electrodes based on distance
% later we have to restrict the neighbouring electrodes to those actually selected in the dataset
channeighbstructmat = (dist<neighbourdist);

% electrode istelf is not a neighbour
channeighbstructmat = (channeighbstructmat .* ~eye(nsensors));

% construct a structured cell array with all neighbours
neighbours=struct;
for i=1:nsensors
  neighbours(i).label       = label{i};
  neighbours(i).neighblabel = label(find(channeighbstructmat(i,:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that computes the neighbourhood geometry from the
% triangulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neighbours] = compneighbstructfromtri(chanpos, label, tri)

nsensors = length(label);

channeighbstructmat = zeros(nsensors,nsensors);
% mark neighbours according to triangulation
for i=1:size(tri, 1)
  channeighbstructmat(tri(i, 1), tri(i, 2)) = 1;
  channeighbstructmat(tri(i, 1), tri(i, 3)) = 1;
  channeighbstructmat(tri(i, 2), tri(i, 1)) = 1;
  channeighbstructmat(tri(i, 3), tri(i, 1)) = 1;
  channeighbstructmat(tri(i, 2), tri(i, 3)) = 1;
  channeighbstructmat(tri(i, 3), tri(i, 2)) = 1;
end

% construct a structured cell array with all neighbours
neighbours = struct;
alldist = [];
for i=1:nsensors
  neighbours(i).label       = label{i};
  neighbidx                 = find(channeighbstructmat(i,:));
  neighbours(i).dist        = sqrt(sum((repmat(chanpos(i, :), numel(neighbidx), 1) - chanpos(neighbidx, :)).^2, 2));
  alldist                   = [alldist; neighbours(i).dist];
  neighbours(i).neighblabel = label(neighbidx);
end

% remove neighbouring channels that are too far away (imporntant e.g. in
% case of missing sensors)
neighbdist = mean(alldist)+3*std(alldist);
for i=1:nsensors
  idx = neighbours(i).dist > neighbdist;
  neighbours(i).dist(idx)         = [];
  neighbours(i).neighblabel(idx)  = [];
end
neighbours = rmfield(neighbours, 'dist');
