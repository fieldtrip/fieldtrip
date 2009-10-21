function neighbours = neighbourselection(cfg,data)

% NEIGHBOURSELECTION finds the neighbours of the channels on the basis of a
% minimum neighbourhood distance (in cfg.neighbourdist). The positions of
% the channel are specified in a gradiometer or electrode configuration or 
% from a layout.
% This configuration can be passed in three ways: 
%  (1) in a configuration field, 
%  (2) in a file whose name is passed in a configuration field, and that can be imported using READ_SENS, or 
%  (3) in a data field. 
%
% Use as
%   neighbours = neighbourselection(cfg, data)
%
% The configuration can contain
%   cfg.neighbourdist = number, maximum distance between neighbouring sensors
%   cfg.elec          = structure with EEG electrode positions
%   cfg.grad          = structure with MEG gradiometer positions
%   cfg.elecfile      = filename containing EEG electrode positions
%   cfg.gradfile      = filename containing MEG gradiometer positions
%   cfg.layout        = filename of the layout, see PREPARE_LAYOUT
%   cfg.feedback      = 'yes' or 'no' (default = 'no')
%
% The following data fields may also be used by NEIGHBOURSELECTION:
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
% $Log: neighbourselection.m,v $
% Revision 1.13  2009/04/07 15:49:40  ingnie
% Only changed help -> fixed typo and mentioned layout possibility in description.
%
% Revision 1.12  2008/11/12 19:22:38  roboos
% documented cfg.layout, added cfg.feedback
%
% Revision 1.11  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.10  2008/03/05 10:46:36  roboos
% moved electrode reading functionality from read_fcdc_elec to read_sens, switched to the use of the new function
%
% Revision 1.9  2008/01/24 17:05:04  roboos
% mention layout in the "undocumented options" section
%
% Revision 1.8  2007/05/14 08:26:31  roboos
% added option to construct neighbours from 2-D layout
%
% Revision 1.7  2006/10/11 09:44:54  roboos
% updated documentation
%
% Revision 1.6  2006/07/12 14:14:59  roboos
% get sens from data.grad/elec
%
% Revision 1.5  2006/07/03 12:57:07  erimar
% Improved help.
%
% Revision 1.4  2006/06/12 08:25:25  erimar
% Added help concerning the structure of cfg.neighbour. Removed
% subfunction involving the calculation of the channel (combination)
% neighourhood geometry.
%
% Revision 1.3  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.2  2006/04/12 09:10:53  roboos
% added a fprintf statement about the number of neighbours that was found
%
% Revision 1.1  2006/04/11 16:15:24  roboos
% created seperate implementation for the construction of the neighbourhood structure, slightly comparable to channelselection
%

fieldtripdefs

% set the defaults
if ~isfield(cfg, 'neighbourdist'),  cfg.neighbourdist = 4; end
if ~isfield(cfg, 'feedback'),       cfg.feedback = 'no';   end

% get the the grad or elec if not present in the data
if isfield(cfg, 'grad')
  fprintf('Obtaining the gradiometer configuration from the configuration.\n');
  sens = cfg.grad;
elseif isfield(cfg, 'elec')
  fprintf('Obtaining the electrode configuration from the configuration.\n');
  sens = cfg.elec;
elseif isfield(cfg, 'gradfile')
  fprintf('Obtaining the gradiometer configuration from a file.\n');
  sens = read_sens(cfg.gradfile);
elseif isfield(cfg, 'elecfile')
  fprintf('Obtaining the electrode configuration from a file.\n');
  sens = read_sens(cfg.elecfile);
elseif isfield(cfg, 'layout')
  fprintf('Using the 2-D layout to determine the neighbours\n');
  lay = prepare_layout(cfg);
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

