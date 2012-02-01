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
%   cfg.method        = 'distance', 'triangulation' or 'template' (default = 'distance')
%   cfg.neighbourdist = number, maximum distance between neighbouring sensors (only for 'distance')
%   cfg.template      = name of the template file, e.g. CTF275_neighb.mat
%   cfg.layout        = filename of the layout, see FT_PREPARE_LAYOUT
%   cfg.elec          = structure with EEG electrode positions
%   cfg.grad          = structure with MEG gradiometer positions
%   cfg.elecfile      = filename containing EEG electrode positions
%   cfg.gradfile      = filename containing MEG gradiometer positions
%   cfg.feedback      = 'yes' or 'no' (default = 'no')
%
% The following data fields may also be used by FT_PREPARE_NEIGHBOURS:
%   data.elec     = structure with EEG electrode positions
%   data.grad     = structure with MEG gradiometer positions
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
% See also FT_NEIGHBOURPLOT, FT_FETCH_SENS

% Copyright (C) 2006-2011, Eric Maris, Jorn M. Horschig, Robert Oostenveld
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', {'method'});

% set the defaults
if ~isfield(cfg, 'feedback'),       cfg.feedback = 'no';         end

hasdata = nargin>1;
if hasdata, data = ft_checkdata(data); end

if strcmp(cfg.method, 'template')
  neighbours = [];
  fprintf('Trying to load sensor neighbours from a template\n');
  if ~isfield(cfg, 'template')
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
  if ~isfield(cfg, 'template')
    if ~isfield(cfg, 'layout')
      error('You need to define a template or layout or give data as an input argument when ft_prepare_neighbours is called with cfg.method=''template''');
    end
    fprintf('Using the 2-D layout filename to determine the template filename\n');
    cfg.template = [strtok(cfg.layout, '.') '_neighb.mat'];
  end
  cfg.template = lower(cfg.template);
  if ~exist(cfg.template, 'file')
    error('Template file could not be found - please check spelling or contact jm.horschig(at)donders.ru.nl if you want to create and share your own template! See also http://fieldtrip.fcdonders.nl/faq/how_can_i_define_my_own_neighbourhood_template');
  end
  load(cfg.template);
  fprintf('Successfully loaded neighbour structure from %s\n', cfg.template);
  
  % only select those channels that are in the data
  if (hasdata)
    neighb_chans = {neighbours(:).label};
    chans = ft_channelselection(data.label, neighb_chans);
    neighb_idx = ismember(neighb_chans, chans);
    neighbours = neighbours(neighb_idx);
  end
else
  % get the the grad or elec if not present in the data
  if hasdata 
    sens = ft_fetch_sens(cfg, data);
  else
    sens = ft_fetch_sens(cfg);
  end
  
  switch lower(cfg.method)
    case 'distance'
      % use a smart default for the distance
      if ~isfield(cfg, 'neighbourdist')
        sens = ft_checkdata(sens, 'hasunits', 'yes');
        if isfield(sens, 'unit') && strcmp(sens.unit, 'm')
          cfg.neighbourdist = 0.04;
        elseif isfield(sens, 'unit') && strcmp(sens.unit, 'dm')
          cfg.neighbourdist = 0.4;
        elseif isfield(sens, 'unit') && strcmp(sens.unit, 'cm')
          cfg.neighbourdist = 4;
        elseif isfield(sens, 'unit') && strcmp(sens.unit, 'mm')
          cfg.neighbourdist = 40;
        else
          % don't provide a default in case the dimensions of the sensor array are unknown
          error('Sensor distance is measured in an unknown unit type');
        end
      end
      
      neighbours = compneighbstructfromgradelec(sens, cfg.neighbourdist);
    case {'triangulation', 'tri'} % the latter for reasons of simplicity
      if size(sens.chanpos, 2)==2 || all(sens.chanpos(:,3)==0)
        % the sensor positions are already on a 2D plane
        prj = sens.chanpos(:,1:2);
      else
        % project sensor on a 2D plane
        prj = elproj(sens.chanpos);
      end
      % make a 2d delaunay triangulation of the projected points
      tri = delaunay(prj(:,1), prj(:,2));
      tri_x = delaunay(prj(:,1)./2, prj(:,2));
      tri_y = delaunay(prj(:,1), prj(:,2)./2);
      tri = [tri; tri_x; tri_y];
      neighbours = compneighbstructfromtri(sens, tri);
    otherwise
      error('Method ''%s'' not known', cfg.method);
  end
end

% removed as from Nov 09 2011 - hope there are no problems with this
% if iscell(neighbours)
%   warning('Neighbourstructure is in old format - converting to structure array');
%   neighbours = fixneighbours(neighbours);
% end

k = 0;
for i=1:length(neighbours)
  k = k + length(neighbours(i).neighblabel);
end
if k==0, error('No neighbours were found!'); end;
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
ft_postamble trackconfig
ft_postamble callinfo
if hasdata
  ft_postamble previous data
end
ft_postamble history neighbours

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that compute the neighbourhood geometry from the
% gradiometer/electrode positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neighbours]=compneighbstructfromgradelec(sens,neighbourdist)

nsensors = length(sens.label);

% compute the distance between all sensors
dist = zeros(nsensors,nsensors);
for i=1:nsensors
  dist(i,:) = sqrt(sum((sens.chanpos(1:nsensors,:) - repmat(sens.chanpos(i,:), nsensors, 1)).^2,2))';
end;

% find the neighbouring electrodes based on distance
% later we have to restrict the neighbouring electrodes to those actually selected in the dataset
channeighbstructmat = (dist<neighbourdist);

% electrode istelf is not a neighbour
channeighbstructmat = (channeighbstructmat .* ~eye(nsensors));

% construct a structured cell array with all neighbours
neighbours=struct;
for i=1:nsensors
  neighbours(i).label       = sens.label{i};
  neighbours(i).neighblabel = sens.label(find(channeighbstructmat(i,:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that computes the neighbourhood geometry from the
% triangulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neighbours]=compneighbstructfromtri(sens,tri)

nsensors = length(sens.label);

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
neighbours=struct;
alldist = [];
for i=1:nsensors
  neighbours(i).label       = sens.label{i};
  neighbidx                 = find(channeighbstructmat(i,:));
  neighbours(i).dist        = sqrt(sum((repmat(sens.chanpos(i, :), numel(neighbidx), 1) - sens.chanpos(neighbidx, :)).^2, 2));
  alldist                   = [alldist; neighbours(i).dist];
  neighbours(i).neighblabel = sens.label(neighbidx);
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
