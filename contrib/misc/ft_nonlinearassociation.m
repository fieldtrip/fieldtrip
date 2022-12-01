function [association] = ft_nonlinearassociation(cfg, data)

% NONLINEARASSOCIATION calculate the association coefficient as a
% function of delay.
%
% In order to estimate the amount of association between all possible
% pairs of MEG sensors the nonlinear association analysis is used.
% It was first developed for EEG data analysis by Pijn and co-workers
% (Lopes da Silva, et al. 1989; Pijn, et al. 1990). The basic principle
% is similar to that of coherence and (cross) correlation, with the
% exception that this nonlinear association method can be applied
% independent of the type of relationship (linear or nonlinear) in
% the data.
% 
% The method is based on the idea that if two signals x and y are
% correlated, a nonlinear regression curve can be calculated that
% represents their relationship. In practice, that regression curve
% is estimated by creating a scatterplot of y versus x, dividing the
% data in segments and describing each segment with a linear regression
% curve. The estimated correlation ratio h2, which gives the reduction
% in variance of y as a result of predicting its values according to
% the regression curve, can be calculated as follows:
% 
% h^2 = (sum(Yi - mean(Y))^2 - sum(Yi - f(Xi))^2) / sum(Yi - mean(Y))^2
% 
% With the sum going over N samples and f(Xi) the estimated value of
% Yi according to the regression line. The h2 coefficient has values
% between 0 (y is completely independent of x) and 1 (y is completely
% determined by x). In the case of a linear relationship between x
% and y, h2 is equal to the well known Pearson correlation coefficient
% (r2). As is the case with cross-correlation, it is possible to
% estimate h2 as a function of time shift () between the signals. The
% h2 is then iteratively calculated for different values of , by
% shifting the signals in comparison to each other, and the value for
% which the maximal h2 is reached can be used as an estimate of the
% time lag between both signals. In deciding what epoch length to use
% in the association analysis, a trade-off has to be made between
% successfully determining the correct delay and h2-value (for which
% large epoch lengths are necessary) and a small enough time-resolution
% (for which small epoch lengths are necessary).
% 
% Use as
%   [association] = ft_nonlinearassociation(cfg, data)
%
% The input data should be organised in a structure as obtained from
% the PREPROCESSING function.
%
% The configuration should contain
%   cfg.channel    = Nx1 cell-array with selection of channels (default = 'all'), see CHANNELSELECTION for details
%   cfg.keeptrials = 'yes' or 'no', process the individual trials or the concatenated data (default = 'no')
%   cfg.trials     = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.fsample    = 1200
%   cfg.maxdelay   = 32/cfg.fsample
%   cfg.delaystep  = 2/cfg.fsample
%   cfg.nr_bins    = 7
%   cfg.offset     = 0
%   cfg.order      = 'Hxy'
%   cfg.timwin     = 0.2
%   cfg.toi        = []
%
% References
% - Lopes da Silva F, Pijn JP, Boeijinga P. (1989): Interdependence of
% EEG signals: linear vs. nonlinear associations and the significance
% of time delays and phase shifts. Brain Topogr 2(1-2):9-18.
% - Pijn JP, Vijn PC, Lopes da Silva FH, Van Ende Boas W, Blanes W.
% (1990): Localization of epileptogenic foci using a new signal
% analytical approach. Neurophysiol Clin 20(1):1-11.

% Copyright (C) 2007, Inge Westmijse

ft_defaults

% set the defaults
if ~isfield(cfg, 'trials'),        cfg.trials  = 'all';             end
if ~isfield(cfg, 'keeptrials'),    cfg.keeptrials  = 'no';          end
if ~isfield(cfg, 'channel'),       cfg.channel = 'all';             end
if ~isfield(cfg, 'toi'),           cfg.toi = [];                    end
if ~isfield(cfg, 'timwin'),        cfg.timwin = 0.2;                end
if ~isfield(cfg, 'fsample'),       cfg.fsample = 1200;              end
if ~isfield(cfg, 'delaystep'),     cfg.delaystep = 2/cfg.fsample;   end
if ~isfield(cfg, 'maxdelay'),      cfg.maxdelay = 32/cfg.fsample;   end
if ~isfield(cfg, 'offset'),        cfg.offset = 0;                  end
if ~isfield(cfg, 'nr_bins'),       cfg.nr_bins = 7;                 end
if ~isfield(cfg, 'order'),         cfg.order = 'Hxy';               end
if ~isfield(cfg, 'feedback'),      cfg.feedback = 'textbar';        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do some bookkeeping on the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
  fprintf('selecting %d trials\n', length(cfg.trials));
  data.trial  = data.trial(cfg.trials);
  data.time   = data.time(cfg.trials);
end
Ntrials  = numel(data.trial);

% select channels of interest
cfg.channel = channelselection(cfg.channel, data.label);
chansel = match_str(data.label, cfg.channel);
fprintf('selecting %d channels\n', length(chansel));
for trial=1:Ntrials
  data.trial{trial} = data.trial{trial}(chansel,:);
end
data.label = data.label(chansel);
Nchans     = length(chansel);

% determine the size of each trial, they can be variable length
Nsamples = zeros(1,Ntrials);
for trial=1:Ntrials
  Nsamples(trial) = size(data.trial{trial},2);
end

if strcmp(cfg.keeptrials, 'no')
  % concatenate all the data into a 2D matrix
  fprintf('concatenating data');
  dat = zeros(Nchans, sum(Nsamples));
  for trial=1:Ntrials
    fprintf('.');
    begsample = sum(Nsamples(1:(trial-1))) + 1;
    endsample = sum(Nsamples(1:trial));
    dat(:,begsample:endsample) = data.trial{trial};
  end
  fprintf('\n');
  fprintf('concatenated data matrix size %dx%d\n', size(dat,1), size(dat,2));
  time = [1/cfg.fsample : size(dat,1)/cfg.fsample : size(dat,1)];
  data.trial = {dat};
  data.time  = {time};
  Ntrials    = 1;
else
  % replace the time axis, since shifted time-axes over trials are not supported
  for i = 1:Ntrials
    data.time{i} = [1/cfg.fsample : 1/cfg.fsample : size(data.trial{i},2) / cfg.fsample];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare all data selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% item 1: channel combinations
numchannelcmb = Nchans*(Nchans-1)/2;
channelcmb = zeros(numchannelcmb, 2);
k = 1;

if (strcmp(cfg.order, 'Hxy'))
  for i=1:Nchans
    for j=(i+1):Nchans
      channelcmb(k, :) = [i j];
      k = k+1;
    end
  end
elseif (strcmp(cfg.order, 'Hyx'))
  for i = Nchans:-1:1
    for j = Nchans : -1:(i+1)
      channelcmb(k,:) = [i j];
      k = k+1;
    end
  end
end
numchannelcmb = size(channelcmb,1);

% item 2: time windows
% this requires
%    cfg.toi = [t1 t2 t3 t4]    center time of each window
%    cfg.timwin = 0.2           in seconds
% TODO implement "time, offset, fsample"

for j = 1:Ntrials
  time = data.time{j};
  numtimwin = size(cfg.toi,2); % initial guess, not all might fit in this trial
  timwin = zeros(numtimwin,2);
  sel = zeros(numtimwin,1);
  for i=1:numtimwin
    timwin(i,1) = cfg.toi(i) - cfg.timwin/2;  % begin of time window
    timwin(i,2) = cfg.toi(i) + cfg.timwin/2;  % end of time window
    sel(i) = timwin(i,1)>=time(1) && timwin(i,2)<=time(end);  % does it fit in?
  end
  timwin        = timwin(find(sel==1),:);
  toi_mat{j}    = cfg.toi(find(sel == 1));                        % update the configuration
  timwin_mat{j} = round((timwin - cfg.offset) * cfg.fsample);     % convert to samples
end

timwin = timwin_mat;
cfg.toi = toi_mat;

% item 3: delays within each timewindow
% this requires
%   cfg.delaystep
%   cfg.maxdelay
delay = -cfg.maxdelay:cfg.delaystep:cfg.maxdelay;
numdelay = length(delay);
% convert to samples
delay = round(delay * cfg.fsample);

if strcmp(cfg.keeptrials, 'yes')
  for trllop=1:Ntrials
    association.trial(trllop) = Do_Association_Calculation( trllop , Ntrials , cfg , timwin , numchannelcmb , numdelay , data , channelcmb , delay );
  end % for each trial
else
  trllop = 1;
  association = Do_Association_Calculation( trllop , Ntrials , cfg , timwin , numchannelcmb , numdelay , data , channelcmb , delay );
end

% add the version details of this function call to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id$';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
association.cfg = cfg;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that does the  computation for each trial
function [association] = Do_Association_Calculation( trllop , Ntrials , cfg , timwin , numchannelcmb , numdelay , data , channelcmb , delay )

fprintf('\nprocessing trial %d from %d\n', trllop, Ntrials);

association = [];
numtimwin = size(timwin{trllop},1);
h2 = zeros(numchannelcmb, numtimwin, numdelay);

l = 0;
maxl = numchannelcmb*numtimwin*numdelay;
progress('init', cfg.feedback, 'Computing nonlinear association for each delay');

for i=1:numchannelcmb
  dat1_chan = data.trial{trllop}(channelcmb(i,1),:);
  dat2_chan = data.trial{trllop}(channelcmb(i,2),:);
  progress(l/maxl,'Computing nonlinear association %d from %d', l, maxl);

  for j=1:numtimwin
    dat1_timwin = dat1_chan(timwin{trllop}(j,1):timwin{trllop}(j,2));
    dat2_timwin = dat2_chan(timwin{trllop}(j,1):timwin{trllop}(j,2));

    for k=1:numdelay
      if delay(k)==0
        % select the complete (unshifted) window
        dat1_delay = dat1_timwin;
        dat2_delay = dat2_timwin;
      elseif delay(k)<0
        % channel 1 is shifted w.r.t. channel 2
        dat1_delay = dat1_timwin((abs(delay(k))+1):end);
        dat2_delay = dat2_timwin(1:(end-abs(delay(k))));
      elseif delay(k)>0
        % channel 2 is shifted w.r.t. channel 1
        dat1_delay = dat1_timwin(1:(end-delay(k)));
        dat2_delay = dat2_timwin((delay(k)+1):end);
      end

      % remove the mean of each snippet of data
      dat1_delay = dat1_delay - mean(dat1_delay);
      dat2_delay = dat2_delay - mean(dat2_delay);

      % do the computation

      [sorted_data1,id] = sort(dat1_delay);
      sorted_data2     = dat2_delay(id);

      data_is = sorted_data1;
      data_js = sorted_data2;

      % divided data_i in bins
      [H,mp_i] = hist(data_is,cfg.nr_bins);

      % Divide data_i and data_j in bins
      bp = 1;
      gem_j = zeros (cfg.nr_bins,1);
      for s = 1:cfg.nr_bins
        gem_j(s) = mean(data_js(bp:bp+H(s)-1));
        bp = bp + H(s);
      end

      % Calculation of line segment and variance per bin
      p=1;
      sp = 1;

      clear unex_var tot_var;

      for u = 1:cfg.nr_bins
        data_is_bin = data_is(sp:sp+H(u)-1)';
        data_js_bin = data_js(sp:sp+H(u)-1)';

        if u == 1
          fx1 = gem_j(u) + (gem_j(u+1) - gem_j(u))/(mp_i(u+1) - mp_i(u))*(data_is_bin(1:round(length(data_is_bin)/2)) - mp_i(u));
        else
          fx1 = gem_j(u-1) + (gem_j(u) - gem_j(u-1))/(mp_i(u) - mp_i(u-1))*(data_is_bin(1:round(length(data_is_bin)/2)) - mp_i(u-1));
        end
        if u == cfg.nr_bins
          fx2 = gem_j(u-1) + (gem_j(u) - gem_j(u-1))/(mp_i(u) - mp_i(u-1))*(data_is_bin(round(length(data_is_bin)/2)+1:end) - mp_i(u-1));
        else
          fx2 = gem_j(u) + (gem_j(u+1) - gem_j(u))/(mp_i(u+1) - mp_i(u))*(data_is_bin(round(length(data_is_bin)/2)+1:end) - mp_i(u));
        end

        % calculation unexplained variance
        ftot = [fx1; fx2];
        unex_var(p:p+length(data_js_bin)-1,:) = (data_js_bin - ftot).^2;

        % calculation total variance
        tot_var(p:p+length(data_js_bin)-1,:) = (data_js_bin - mean(data_js)).^2;

        p = p+length(data_is_bin);
        sp = sp + H(u);
      end

      % Calculation of association coefficient h2
      h2(i,j,k) = ((sum(tot_var) - sum(unex_var)) / sum(tot_var))*100;
      l=l+1;
    end
  end
end

progress('close');

% collect the results
if (size(h2, 1) ~= numchannelcmb)
  error('number of channel combinations is not right');
end
if (size(h2, 2) ~= numtimwin)
  error('number of timewindows does not match input');
end
if (size(h2, 3) ~= numdelay)
  error('number of delays does not match input');
end
h2 = reshape(h2, numchannelcmb*numtimwin, numdelay);
[m, indx] = max(h2, [], 2);
h2 = reshape(m, numchannelcmb, numtimwin);
d = (delay(indx)./cfg.fsample)*1000; % convert to miliseconds (to compare with original method)
d = reshape(d, numchannelcmb, numtimwin);

% FIXME this does not work: one is not supposed to rely on data.cfg.trl, use data.time please!
% association.time = cfg.toi{trllop} + data.cfg.trl(trllop,1);

association.h2 = h2;
association.delay = d;

return % from Do_Association_Calculation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that computes the mean over all values without checking the dimensions
% this function is approximately 2x faster on 1 dimensional data than the matlab version
function y = mean(x)
y = sum(x(:))./numel(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that bins the elements of Y into equally spaced containers
% this function is approximately 2x faster on 1 dimensional data than the matlab version
function [nn, x] = hist(y, x)
miny = min(y);
maxy = max(y);
binwidth = (maxy - miny) ./ x;
xx = miny:binwidth:maxy;
x  = xx(1:(end-1)) + binwidth/2;
% Shift bins so the interval is ( ] instead of [ ).
xx = full(real(xx)); y = full(real(y)); % For compatibility
bins = xx + eps(xx);
nn = histc(y',[-inf bins],1);
% Combine first bin with 2nd bin and last bin with next to last bin
nn(2,:) = nn(2,:)+nn(1,:);
nn(end-1,:) = nn(end-1,:)+nn(end,:);
nn = nn(2:end-1,:);

