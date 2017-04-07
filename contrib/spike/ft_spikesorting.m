function [spike] = ft_spikesorting(cfg, spike)

% FT_SPIKESORTING performs clustering of spike-waveforms and returns the
% unit number to which each spike belongs.
%
% Use as
%   [spike] = ft_spikesorting(cfg, spike)
%
% The configuration can contain
%   cfg.channel         cell-array with channel selection (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.method          'kmeans', 'ward'
%   cfg.feedback        'no', 'text', 'textbar', 'gui' (default = 'textbar')
%   cfg.kmeans          substructure with additional low-level options for this method
%   cfg.ward            substructure with additional low-level options for this method
%   cfg.ward.distance   'L1', 'L2', 'correlation', 'cosine'
%
% The input spike structure can be imported using READ_FCDC_SPIKE and should contain
%   spike.label     = 1 x Nchans cell-array, with channel labels
%   spike.waveform  = 1 x Nchans cell-array, each element contains a matrix (Nsamples x Nspikes), can be empty
%   spike.timestamp = 1 x Nchans cell-array, each element contains a vector (1 x Nspikes)
%   spike.unit      = 1 x Nchans cell-array, each element contains a vector (1 x Nspikes)
%
% See also FT_READ_SPIKE, FT_SPIKEDOWNSAMPLE

% Copyright (C) 2006-2007, Robert Oostenveld
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
ft_preamble provenance spike
ft_preamble trackconfig

% these were supported in the past, but are not any more (for consistency with other spike functions)
cfg = ft_checkconfig(cfg, 'forbidden', 'inputfile', ...
                                       'outputfile');  % see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1056

% set the defaults
if ~isfield(cfg, 'feedback'),       cfg.feedback = 'textbar';    end
if ~isfield(cfg, 'method'),         cfg.method = 'ward';         end
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';         end

if isequal(cfg.method, 'ward')
  if ~isfield(cfg, 'ward'),           cfg.ward          = [];       end
  if ~isfield(cfg.ward, 'distance'),  cfg.ward.distance = 'L2';     end
  if ~isfield(cfg.ward, 'optie'),     cfg.ward.optie    = 'ward';   end
  if ~isfield(cfg.ward, 'aantal'),    cfg.ward.aantal   = 10;       end
  if ~isfield(cfg.ward, 'absoluut'),  cfg.ward.absoluut = 1;        end
end

if isequal(cfg.method, 'kmeans')
  if ~isfield(cfg, 'kmeans'),         cfg.kmeans        = [];       end
  if ~isfield(cfg.kmeans, 'aantal'),  cfg.kmeans.aantal = 10;       end
end

% select the channels
cfg.channel = ft_channelselection(cfg.channel, spike.label);
sel = match_str(spike.label, cfg.channel);
spike.label     = spike.label(sel);
spike.waveform  = spike.waveform(sel);
spike.timestamp = spike.timestamp(sel);
spike.unit      = {};  % this will be assigned by the clustering

nchan = length(spike.label);
for chanlop=1:nchan
  label     = spike.label{chanlop};
  waveform  = spike.waveform{chanlop};
  timestamp = spike.timestamp{chanlop};
  nspike    = size(waveform,2);
  unit      = zeros(1,nspike);
  fprintf('sorting %d spikes in channel %s\n', nspike, label);
  
  switch cfg.method
    case 'ward'
      dist = ward_distance(cfg, waveform);
      [grootte,ordening] = cluster_ward(dist,cfg.ward.optie,cfg.ward.absoluut,cfg.ward.aantal);
      for i=1:cfg.ward.aantal
        begsel = sum([1 grootte(1:(i-1))]);
        endsel = sum([0 grootte(1:(i  ))]);
        unit(ordening(begsel:endsel)) = i;
      end
      
    case 'kmeans'
      unit = kmeans(waveform', cfg.kmeans.aantal)';
      
      %               'sqEuclidean'  - Squared Euclidean distance
      %               'cityblock'    - Sum of absolute differences, a.k.a. L1
      %               'cosine'       - One minus the cosine of the included angle
      %                                between points (treated as vectors)
      %               'correlation'  - One minus the sample correlation between
      %                                points (treated as sequences of values)
      
    otherwise
      error('unsupported clustering method');
  end
  
  % remember the sorted units
  spike.unit{chanlop} = unit;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous   spike
ft_postamble provenance spike
ft_postamble history    spike

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that computes the distance between all spike waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist] = ward_distance(cfg, waveform)
nspike = size(waveform,2);
dist = zeros(nspike, nspike);
ft_progress('init', cfg.feedback, 'computing distance');
switch lower(cfg.ward.distance)
  case 'l1'
    for i=1:nspike
      ft_progress((i-1)/nspike, 'computing distance for spike %d/%d', i, nspike);
      for j=2:nspike
        dist(i,j) = sum(abs(waveform(:,j)-waveform(:,i)));
        dist(j,i) = dist(i,j);
      end
    end
  case 'l2'
    for i=1:nspike
      ft_progress((i-1)/nspike, 'computing distance for spike %d/%d', i, nspike);
      for j=2:nspike
        dist(i,j) = sqrt(sum((waveform(:,j)-waveform(:,i)).^2));
        dist(j,i) = dist(i,j);
      end
    end
  case 'correlation'
    for i=1:nspike
      ft_progress((i-1)/nspike, 'computing distance for spike %d/%d', i, nspike);
      for j=2:nspike
        dist(i,j) = corrcoef(waveform(:,j),waveform(:,i));
        dist(j,i) = dist(i,j);
      end
    end
  case 'cosine'
    for i=1:nspike
      ft_progress((i-1)/nspike, 'computing distance for spike %d/%d', i, nspike);
      for j=2:nspike
        x = waveform(:,j);
        y = waveform(:,i);
        x = x./norm(x);
        y = y./norm(y);
        dist(i,j) = dot(x, y);
        dist(j,i) = dist(i,j);
      end
    end
    
  otherwise
    error('unsupported distance metric');
end
ft_progress('close');
