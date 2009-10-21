function [spike] = spikesorting(cfg, spike);

% SPIKESORTING performs clustering of spike-waveforms and returns the unit number to which each spike belongs.
%
% Use as
%   [spike] = spikesorting(cfg, spike)
%
% The configuration can contain
%   cfg.channel         cell-array with channel selection (default = 'all'), see CHANNELSELECTION for details
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
% See also READ_FCDC_SPIKE, SPIKEDOWNSAMPLE

% Copyright (C) 2006-2007, Robert Oostenveld
%
% $Log: spikesorting.m,v $
% Revision 1.5  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.4  2007/05/02 16:04:30  roboos
% added implementation of kmeans clustering, changed the function to reflect the updated data type
%
% Revision 1.3  2007/01/09 09:53:31  roboos
% removed neuralynx_ds, added plexon_plx
% the code currently is broken due to the changed behaviour of the reading functions
%
% Revision 1.2  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.1  2006/03/10 13:43:08  roboos
% initial implementation, only ward clustering sofar
%

fieldtripdefs

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
cfg.channel = channelselection(cfg.channel, spike.label);
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

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: spikesorting.m,v 1.5 2008/09/22 20:17:44 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous    = spike.cfg;     end
% remember the configuration
spike.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that computes the distance between all spike waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist] = ward_distance(cfg, waveform);
nspike = size(waveform,2);
dist = zeros(nspike, nspike);
progress('init', cfg.feedback, 'computing distance');
switch lower(cfg.ward.distance)
  case 'l1'
    for i=1:nspike
      progress((i-1)/nspike, 'computing distance for spike %d/%d', i, nspike);
      for j=2:nspike
        dist(i,j) = sum(abs(waveform(:,j)-waveform(:,i)));
        dist(j,i) = dist(i,j);
      end
    end
  case 'l2'
    for i=1:nspike
      progress((i-1)/nspike, 'computing distance for spike %d/%d', i, nspike);
      for j=2:nspike
        dist(i,j) = sqrt(sum((waveform(:,j)-waveform(:,i)).^2));
        dist(j,i) = dist(i,j);
      end
    end
  case 'correlation'
    for i=1:nspike
      progress((i-1)/nspike, 'computing distance for spike %d/%d', i, nspike);
      for j=2:nspike
        dist(i,j) = corrcoef(waveform(:,j),waveform(:,i));
        dist(j,i) = dist(i,j);
      end
    end
  case 'cosine'
    for i=1:nspike
      progress((i-1)/nspike, 'computing distance for spike %d/%d', i, nspike);
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
progress('close');
