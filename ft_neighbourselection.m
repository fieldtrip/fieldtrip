function neighbours = ft_neighbourselection(cfg,data)

% FT_NEIGHBOURSELECTION finds the neighbours of the channels on the basis of a
% minimum neighbourhood distance (in cfg.neighbourdist). The positions of
% the channel are specified in a gradiometer or electrode configuration or
% from a layout.
% This configuration can be passed in three ways:
%  (1) in a configuration field,
%  (2) in a file whose name is passed in a configuration field, and that can be imported using READ_SENS, or
%  (3) in a data field.
%
% Use as
%   neighbours = ft_neighbourselection(cfg, data)
%
% The configuration can contain
%   cfg.neighbourdist = number, maximum distance between neighbouring sensors
%   cfg.elec          = structure with EEG electrode positions
%   cfg.grad          = structure with MEG gradiometer positions
%   cfg.elecfile      = filename containing EEG electrode positions
%   cfg.gradfile      = filename containing MEG gradiometer positions
%   cfg.layout        = filename of the layout, see FT_PREPARE_LAYOUT
%   cfg.feedback      = 'yes' or 'no' (default = 'no')
%
% The following data fields may also be used by FT_NEIGHBOURSELECTION:
%   data.elec     = structure with EEG electrode positions
%   data.grad     = structure with MEG gradiometer positions
%
% The output:
%   neighbours     = definition of neighbours for each channel,
%     which is structured like this:
%        neighbours{1}.label = 'Fz';
%        neighbours{1}.neighblabel = {'Cz', 'F3', 'F3A', 'FzA', 'F4A', 'F4'};
%        neighbours{2}.label = 'Cz';
%        neighbours{2}.neighblabel = {'Fz', 'F4', 'RT', 'RTP', 'P4', 'Pz', 'P3', 'LTP', 'LT', 'F3'};
%        neighbours{3}.label = 'Pz';
%        neighbours{3}.neighblabel = {'Cz', 'P4', 'P4P', 'Oz', 'P3P', 'P3'};
%        etc.
%        (Note that a channel is not considered to be a neighbour of itself.)

% Copyright (C) 2006-2008, Eric Maris, Robert Oostenveld
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

ft_defaults

% set the defaults
if ~isfield(cfg, 'feedback'),       cfg.feedback = 'no';   end

% the default value for this is set later
% if ~isfield(cfg, 'neighbourdist'),  cfg.neighbourdist = 4; end

% get the the grad or elec if not present in the data
if isfield(cfg, 'grad')
  fprintf('Obtaining the gradiometer configuration from the configuration.\n');
  sens = cfg.grad;
elseif isfield(cfg, 'elec')
  fprintf('Obtaining the electrode configuration from the configuration.\n');
  sens = cfg.elec;
elseif isfield(cfg, 'gradfile')
  fprintf('Obtaining the gradiometer configuration from a file.\n');
  sens = ft_read_sens(cfg.gradfile);
elseif isfield(cfg, 'elecfile')
  fprintf('Obtaining the electrode configuration from a file.\n');
  sens = ft_read_sens(cfg.elecfile);
elseif isfield(cfg, 'layout')
  fprintf('Using the 2-D layout to determine the neighbours\n');
  lay = ft_prepare_layout(cfg);
  sens = [];
  sens.label = lay.label;
  sens.pnt = lay.pos;
  sens.pnt(:,3) = 0;
elseif isfield(data, 'grad')
  fprintf('Using the gradiometer configuration from the dataset.\n');
  sens = data.grad;
elseif isfield(data, 'elec')
  fprintf('Using the electrode configuration from the dataset.\n');
  sens = data.elec;
end
if ~isstruct(sens)
  error('Did not find gradiometer or electrode information.');
end;

% use a smart default for the distance
if ~isfield(cfg, 'neighbourdist'),
  if     isfield(sens, 'unit') && strcmp(sens.unit, 'cm')
    cfg.neighbourdist = 4;
  elseif isfield(sens, 'unit') && strcmp(sens.unit, 'mm')
    cfg.neighbourdist = 40;
  else
    % don't provide a default in case the dimensions of the sensor array are unknown
  end
end

neighbours = compneighbstructfromgradelec(sens, cfg.neighbourdist);

k = 0;
for i=1:length(neighbours)
  k = k + length(neighbours{i}.neighblabel);
end
fprintf('there are on average %.1f neighbours per channel\n', k/length(neighbours));

if strcmp(cfg.feedback, 'yes')
  % give some graphical feedback
  if all(sens.pnt(:,3)==0)
    % the sensor positions are already projected on a 2D plane
    proj = sens.pnt(:,1:2);
  else
    % project the 3D positions onto a 2D plane
    proj = elproj(sens.pnt);
  end
  figure
  for i=1:length(neighbours)
    cla
    this = neighbours{i};
    sel1 = match_str(sens.label, this.label);
    sel2 = match_str(sens.label, this.neighblabel);
    plot(proj(:,1), proj(:,2), 'k.');
    axis equal
    title(this.label);
    axis off
    for j=1:length(this.neighblabel)
      x1 = proj(sel1,1);
      y1 = proj(sel1,2);
      x2 = proj(sel2(j),1);
      y2 = proj(sel2(j),2);
      X = [x1 x2];
      Y = [y1 y2];
      line(X, Y, 'color', 'r');
    end
    drawnow
    pause(0.1);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that compute the neighbourhood geometry from the
% gradiometer/electrode positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neighbours]=compneighbstructfromgradelec(sens,neighbourdist)

nsensors = length(sens.label);

% compute the distance between all sensors
dist = zeros(nsensors,nsensors);
for i=1:nsensors
  dist(i,:) = sqrt(sum((sens.pnt(1:nsensors,:) - repmat(sens.pnt(i,:), nsensors, 1)).^2,2))';
end;

% find the neighbouring electrodes based on distance
% later we have to restrict the neighbouring electrodes to those actually selected in the dataset
channeighbstructmat = (dist<neighbourdist);

% electrode istelf is not a neighbour
channeighbstructmat = (channeighbstructmat .* ~eye(nsensors));

% construct a structured cell array with all neighbours
neighbours=cell(1,nsensors);
for i=1:nsensors
  neighbours{i}.label       = sens.label{i};
  neighbours{i}.neighblabel = sens.label(find(channeighbstructmat(i,:)));
end

