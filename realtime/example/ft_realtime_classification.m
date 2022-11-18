function ft_realtime_classification(cfg)

% FT_REALTIME_CLASSIFICATION is an example realtime application for online
% classification of the data. It should work both for EEG and MEG.
%
% Use as
%   ft_realtime_classification(cfg)
% with the following configuration options
%   cfg.channel    = cell-array, see FT_CHANNELSELECTION (default = 'all')
%   cfg.trialfun   = string with the trial function
%
% The source of the data is configured as
%   cfg.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.datafile      = string
%   cfg.headerfile    = string
%   cfg.eventfile     = string
%   cfg.dataformat    = string, default is determined automatic
%   cfg.headerformat  = string, default is determined automatic
%   cfg.eventformat   = string, default is determined automatic
%
% This function works with two-class data that is timelocked to a trigger.
% Data selection is based on events that should be present in the
% datastream or datafile. The user should specify a trial function that
% selects pieces of data to be classified, or pieces of data on which the
% classifier has to be trained.The trialfun should return segments in a
% trial definition (see FT_DEFINETRIAL). The 4th column of the trl matrix
% should contain the class label (number 1 or 2). The 5th colum of the trl
% matrix should contain a flag indicating whether it belongs to the test or
% to the training set (0 or 1 respectively).
%
% Example usage:
%   cfg = [];
%   cfg.dataset  = 'Subject01.ds';
%   cfg.trialfun = 'trialfun_Subject01';
%   ft_realtime_classification(cfg);
%
% To stop the realtime function, you have to press Ctrl-C

% Undocumented options:
%   cfg.timeout = scalar, time in seconds after which the function stops.
%                 Default value is inf, but may be set to a finite number
%                 (so that it stops executing when running without user
%                 interaction).

% Copyright (C) 2009, Robert Oostenveld
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

% this makes use of an external classification toolbox
ft_hastoolbox('prtools', 1);

% set the default configuration options
cfg.dataformat   = ft_getopt(cfg, 'dataformat',   []); % default is detected automatically
cfg.headerformat = ft_getopt(cfg, 'headerformat', []); % default is detected automatically
cfg.eventformat  = ft_getopt(cfg, 'eventformat', []);  % default is detected automatically
cfg.channel      = ft_getopt(cfg, 'channel',    'all');
cfg.bufferdata   = ft_getopt(cfg, 'bufferdata', 'last'); % first or last
cfg.timeout      = ft_getopt(cfg, 'timeout',    inf);

% translate dataset into datafile+headerfile
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
cfg = ft_checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear ft_read_header
% start by reading the header from the realtime buffer
hdr = ft_read_header(cfg.headerfile, 'cache', true);

% define a subset of channels for reading
cfg.channel = ft_channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);
nchan       = length(chanindx);

if nchan==0
  ft_error('no channels were selected');
end

% these are for the data handling
prevSample = 0;
count      = 0;

% measure the timeing
tic;
t(1) = toc;
s(1) = 0;

% these are for the classification
W           = [];
correct     = [];
train_class = [];
train_dat   = [];
clear(cfg.trialfun);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while t(end)<cfg.timeout

  % determine latest header and event information
  event     = ft_read_event(cfg.dataset, 'minsample', prevSample+1);  % only consider events that are later than the data processed sofar
  hdr       = ft_read_header(cfg.dataset, 'cache', true);             % the trialfun might want to use this, but it is not required
  cfg.event = event;                                               % store it in the configuration, so that it can be passed on to the trialfun
  cfg.hdr   = hdr;                                                 % store it in the configuration, so that it can be passed on to the trialfun

  % evaluate the trialfun, note that the trialfun should not re-read the events and header
  fprintf('evaluating ''%s'' based on %d events\n', cfg.trialfun, length(event));
  trl = feval(cfg.trialfun, cfg);

  % the code below assumes that the 4th column of the trl matrix contains
  % the class label and the 5th column a boolean indicating whether it is a
  % training set item or test set item
  if size(trl,2)<4
    trl(:,4) = nan; % don't asign a default class
  end
  if size(trl,2)<5
    trl(:,5) = 0; % assume that it is a test set item
  end

  fprintf('processing %d trials\n', size(trl,1));

  for trllop=1:size(trl,1)

    begsample = trl(trllop,1);
    endsample = trl(trllop,2);
    class     = trl(trllop,4);
    train     = trl(trllop,5)==1;
    test      = trl(trllop,5)==0;

    % remember up to where the data was read
    prevSample  = endsample;
    count       = count + 1;
    fprintf('-------------------------------------------------------------------------------------\n');
    fprintf('processing segment %d from sample %d to %d, class = %d, train = %d\n', count, begsample, endsample, class, train);

    % read data segment from buffer
    dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % keep track of the timing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t(end+1) = toc;
    s(end+1) = endsample;

    % compute the cummulative and instantaneous number of samples per second
    % compare these to the sampling frequency to get the relative acceleration factor
    instantaneous = [nan diff(s) ./ diff(t)];
    cumulative    = (s-s(1)) ./ (t-t(1));
    semilogy([instantaneous(:) cumulative(:)]/hdr.Fs, '.');
    title('acceleration factor');
    legend({'instantaneous', 'cumulative'});
    % force Matlab to update the figure
    drawnow

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to the processing of the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % apply some preprocessing options
    dat = ft_preproc_baselinecorrect(dat);

    if test
      % retrain the classifier based on the accumulated training data
      if isempty(W) && numel(unique(train_class))==2
        % only if the classifier needs to be retrained and can be retrained
        fprintf('retraining the classifier based on %d examples\n', length(train_class));
        A = dataset(train_dat, train_class);
        W = svc(A);
      end

      % classify this trial
      if ~isempty(W)
        [nchan, nsmp] = size(dat);
        dat = reshape(dat, [1, nchan*nsmp]);
        B   = dataset(dat, class);
        Bc  = B*W;
        estimate = labeld(Bc);          % this is the estimated class
      else
        ft_warning('classifier has not yet been trained');
        estimate = nan;
      end

      % keep track of the classification performance
      fprintf('estimated class = %d, real class = %d\n', estimate, class);
      if ~isnan(class)
        % this can only be done if the true class is known
        correct(end+1) = (estimate==class);
        fprintf('classification rate = %d%%\n', round(mean(correct)*100));
      end
    end % if test

    if train
      % delete the previously trained classifier
      W = [];
      % add the current trial to the training data
      fprintf('adding one example to the training dataset\n');
      [nchan, nsmp] = size(dat);
      dat = reshape(dat, [1, nchan*nsmp]);
      if isempty(train_dat)
        train_dat   = dat;
        train_class = class;
      else
        train_dat   = cat(1, train_dat,   dat);
        train_class = cat(1, train_class, class);
      end
    end % if train
  end % looping over new trials
end % while true

