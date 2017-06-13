function [cfg, spike] = ft_spikedetection(cfg)

% FT_SPIKEDETECTION reads continuous data from disk and detects spikes. The
% function writes the unsorted spike waveforms to disk in another file.
%
% Use as
%   cfg = ft_spikedetection(cfg)
%
% The configuration options can contain
%   cfg.dataset             = string with the input dataset
%   cfg.output              = string with the output dataset (default is determined automatic)
%   cfg.dataformat          = string with the output dataset format, see FT_WRITE_FCDC_SPIKE
%   cfg.method              = string with the method to use, can be 'all', 'zthresh', 'ztrig', 'flank'
%   cfg.interactive         = 'yes' or 'no'
%   cfg.timestampdefinition = 'orig' or 'sample'
%
% The default is to process the full dataset. You can select a latency range with
%   cfg.latency          = [begin end], default is [0 inf]
% or you can specify multiple latency segments with
%   cfg.latency          = [b1 e1; b2 e2; ...]
%
% Specific settings for the zthresh spike detection method are
%   cfg.zthresh.neg      = negative threshold, e.g. -3
%   cfg.zthresh.pos      = positive threshold, e.g.  3
%   cfg.zthresh.offset   = number of samples before peak (default = 16)
%   cfg.zthresh.mindist  = mininum distance in samples between  detected peaks
%
% Specific settings for the flank spike detection method are
%   cfg.flank.value      = positive or negative threshold
%   cfg.flank.offset     = number of samples before peak
%   cfg.flank.ztransform = 'yes' or 'no'
%   cfg.flank.mindist    = mininum distance in samples between  detected peaks
%
% Furthermore, the configuration can contain options for preprocessing
%   cfg.preproc.lpfilter      = 'no' or 'yes'  lowpass filter
%   cfg.preproc.hpfilter      = 'no' or 'yes'  highpass filter
%   cfg.preproc.bpfilter      = 'no' or 'yes'  bandpass filter
%   cfg.preproc.lnfilter      = 'no' or 'yes'  line noise removal using notch filter
%   cfg.preproc.dftfilter     = 'no' or 'yes'  line noise removal using discrete fourier transform
%   cfg.preproc.medianfilter  = 'no' or 'yes'  jump preserving median filter
%   cfg.preproc.lpfreq        = lowpass  frequency in Hz
%   cfg.preproc.hpfreq        = highpass frequency in Hz
%   cfg.preproc.bpfreq        = bandpass frequency range, specified as [low high] in Hz
%   cfg.preproc.lnfreq        = line noise frequency in Hz, default 50Hz
%   cfg.preproc.lpfiltord     = lowpass  filter order
%   cfg.preproc.hpfiltord     = highpass filter order
%   cfg.preproc.bpfiltord     = bandpass filter order
%   cfg.preproc.lnfiltord     = line noise notch filter order
%   cfg.preproc.medianfiltord = length of median filter
%   cfg.preproc.lpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.preproc.hpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.preproc.bpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.preproc.lpfiltdir     = filter direction, 'twopass' (default) or 'onepass'
%   cfg.preproc.hpfiltdir     = filter direction, 'twopass' (default) or 'onepass'
%   cfg.preproc.bpfiltdir     = filter direction, 'twopass' (default) or 'onepass'
%   cfg.preproc.detrend       = 'no' or 'yes'
%   cfg.preproc.demean        = 'no' or 'yes'
%   cfg.preproc.baselinewindow = [begin end] in seconds, the default is the complete trial
%   cfg.preproc.hilbert       = 'no' or 'yes'
%   cfg.preproc.rectify       = 'no' or 'yes'

% Copyright (C) 2005-2008, Robert Oostenveld
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
ft_preamble provenance
ft_preamble trackconfig

% set the general defaults
if ~isfield(cfg, 'dataset'),            cfg.dataset = [];             end
if ~isfield(cfg, 'output'),             cfg.output = [];              end
if ~isfield(cfg, 'channel'),            cfg.channel = 'all';          end
if ~isfield(cfg, 'channelprefix'),      cfg.channelprefix = [];       end
if ~isfield(cfg, 'latency'),            cfg.latency = [0 inf];        end
if ~isfield(cfg, 'dataformat'),         cfg.dataformat = [];          end
if ~isfield(cfg, 'headerformat'),       cfg.headerformat = [];        end
% set the specific defaults
if ~isfield(cfg, 'method'),             cfg.method = 'zthresh';       end
if ~isfield(cfg, 'adjustselection'),    cfg.adjustselection = 'yes';  end
if ~isfield(cfg, 'interactive'),        cfg.interactive = 'no';       end
if ~isfield(cfg, 'chanvals'),           cfg.chanvals = [];            end

% set the defaults for the various spike detection methods
switch cfg.method
  case 'all'
    if ~isfield(cfg, 'all'),              cfg.all = [];                  end
  case 'zthresh'
    if ~isfield(cfg, 'zthresh'),          cfg.zthresh = [];              end
    if ~isfield(cfg.zthresh, 'neg'),      cfg.zthresh.neg     = -3;      end
    if ~isfield(cfg.zthresh, 'pos'),      cfg.zthresh.pos     =  3;      end
    if ~isfield(cfg.zthresh, 'offset'),   cfg.zthresh.offset  = -16;     end % in samples
    if ~isfield(cfg.zthresh, 'mindist'),  cfg.zthresh.mindist =  0;      end % in samples
  case 'flank'
    if ~isfield(cfg, 'flank'),            cfg.flank = [];                end
    if ~isfield(cfg.flank, 'ztransform'), cfg.flank.ztransform = 'yes';  end
    if ~isfield(cfg.flank, 'value'),      cfg.flank.value = 1.5;         end % trigger threshold value
    if ~isfield(cfg.flank, 'offset'),     cfg.flank.offset = 6;          end % in samples
    if ~isfield(cfg.flank, 'mindist'),    cfg.flank.mindist = 0;         end % in samples
  otherwise
    error('unsupported option for cfg.method');
end

% ensure that the preproc specific options are located in the preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});
cfg.preproc = ft_checkconfig(cfg.preproc, 'renamed', {'blc', 'demean'});
cfg.preproc = ft_checkconfig(cfg.preproc, 'renamed', {'blcwindow', 'baselinewindow'});

status = mkdir(cfg.output);
if ~status
  error('error creating spike output dataset %s', cfg.output);
end

% read the header of the completete dataset
hdr = ft_read_header(cfg.dataset, 'headerformat', cfg.headerformat);
cfg.channel = ft_channelselection(cfg.channel, hdr.label);
chansel = match_str(hdr.label, cfg.channel);

if strcmp(cfg.timestampdefinition, 'sample')
  % the default would be to keep the original definition of timestamps as determined from looking at the file
  % here the definition of timestamps is changed to correspond with samples at the original sampling rate
  hdr.TimeStampPerSample = 1;
  hdr.FirstTimeStamp     = 1;
  hdr.LastTimeStamp      = hdr.nSamples*hdr.nTrials;
end

if hdr.nSamples<1
  error('the input dataset contains no samples');
elseif length(chansel)<1
  error('the input selection contains no channels');
end

% give some feedback, based on the complete data
fprintf('data contains %10d channels\n', hdr.nChans);
fprintf('selected      %10d channels\n', length(chansel));
numsample = [];
numsegment = size(cfg.latency,1);
for j=1:numsegment
  begsample(j) = max(round(cfg.latency(j,1) * hdr.Fs + 1), 1);
  endsample(j) = min(round(cfg.latency(j,2) * hdr.Fs    ), hdr.nSamples);
  numsample(j) = endsample(j) - begsample(j) + 1;
  cfg.latency(j,1) = (begsample(j)-1)/hdr.Fs;
  cfg.latency(j,2) = (endsample(j)  )/hdr.Fs;
end
numsample = sum(numsample);
fprintf('data contains %10d samples\n', hdr.nSamples);
fprintf('selected      %10d samples in %d segments\n', numsample, numsegment);

s = floor(hdr.nSamples ./ hdr.Fs);
m = floor(s/60);
h = floor(m/60);
m = m - 60*h;
s = s - 60*m - 60*60*h;
fprintf('duration of data      %02dh:%02dm:%02ds\n', h, m, s);

s = floor(numsample ./ hdr.Fs);
m = floor(s/60);
h = floor(m/60);
m = m - 60*h;
s = s - 60*m - 60*60*h;
fprintf('duration of selection %02dh:%02dm:%02ds\n', h, m, s);

fprintf('estimated memory usage %d MB\n', round((numsample*(8+8+2))/(1024^2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process each channel separetely
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=chansel(:)'
  fprintf('processing channel %d, ''%s''\n', i, hdr.label{i});

  % remain in the interactive phase as long as the user desires
  % the interactive loop is also used once when imediately writing to file
  runloop = true;
  newdata = true;

  while runloop
    % the loop is used once when writing to file, or multiple times for interactive use
    runloop = false;

    % reading and filtering may be repeatedly done if the user interactively specified another latency selection
    if newdata
      fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
      numsegment = size(cfg.latency,1);
      clear begsample endsample numsample
      for j=1:numsegment
        begsample(j) = max(round(cfg.latency(j,1) * hdr.Fs + 1), 1);
        endsample(j) = min(round(cfg.latency(j,2) * hdr.Fs    ), hdr.nSamples);
        numsample(j) = endsample(j) - begsample(j) + 1;
        cfg.latency(j,1) = (begsample(j)-1)/hdr.Fs;
        cfg.latency(j,2) = (endsample(j)  )/hdr.Fs;
      end
      % read from a single channel and concatenate the segments into one vector
      org = zeros(1,sum(numsample));
      for j=1:numsegment
        fprintf('reading channel %s, latency from %f to %f\n', hdr.label{i}, cfg.latency(j,1), cfg.latency(j,2));
        buf = ft_read_data(cfg.dataset, 'header', hdr, 'begsample', begsample(j), 'endsample', endsample(j), 'chanindx', i, 'dataformat', cfg.dataformat);
        if j==1
          begsegment = 1;
          endsegment = numsample(j);
        else
          begsegment = sum(numsample(1:j-1)) + 1;
          endsegment = sum(numsample(1:j  ))    ;
        end
        % concatenate the data into one large vector
        org(begsegment:endsegment) = buf;
        clear buf;
      end
      % apply preprocessing
      fprintf('applying preprocessing options\n');
      dat = preproc(org, hdr.label(i), offset2time(0, hdr.Fs, size(org,2)), cfg.preproc);      
    end % if newdata

    peaks = [];

    % this while loop is used to automatically adjust the spike threshold
    % if too few/ too many waveforms were selected
    % it will break out of the while loop at the end
    numadjustment = 1;
    while (1)

      switch cfg.method
        case 'all'
          % place a spike marker at each 32 samples in the signal
          % this makes a lot of waveforms that together reconstruct the complete continuous data
          for j=1:32:(length(dat)-32)
            peaks = [peaks j];
          end

        case 'zthresh'
          % do peak detection on z-transformed signal
          zdat = (dat-mean(dat))./std(dat);
          if ~isinf(cfg.zthresh.neg)
            % mindist is expressed in samples rather than sec
            dum = peakdetect3(-zdat, -cfg.zthresh.neg, cfg.zthresh.mindist);
            peaks = [peaks dum];
          end
          if ~isinf(cfg.zthresh.pos)
            % mindist is expressed in samples rather than sec
            dum   = peakdetect3( zdat,  cfg.zthresh.pos, cfg.zthresh.mindist);
            peaks = [peaks dum];
          end
          % the mindist is not honored for spikes followed immediately by a spike of the other sign
          peaks = sort(peaks);
          % the begin of each waveform is shifted, to ensure that the spike is in the middle
          peaks = peaks - cfg.zthresh.offset;
          % ensure that no spikes are found within the first and last segment of the data
          % since it is not possible to extract complete waveforms there
          peaks(find(peaks<32))               = [];
          peaks(find((length(dat)-peaks)<32)) = [];

          % --- store the thres val for current channnel
          cfg.chanvals(find(chansel==i),1:3) = [cfg.zthresh.neg cfg.zthresh.pos cfg.zthresh.mindist];

        case 'flank'
          if strcmp(cfg.flank.ztransform, 'yes')
            zdat = (dat-mean(dat))./std(dat);
          else
            zdat = dat;
          end
          if cfg.flank.value>0
            peaks = find(diff(zdat>cfg.flank.value)==1);
          elseif cfg.flank.value<0
            peaks = find(diff(zdat<cfg.flank.value)==-1);
          end

          % prevent thres crossing within mindist samples
          if ~isempty(peaks) && ~isempty(cfg.flank.mindist)
            pd = [inf diff(peaks)];
            peaks = peaks(pd>cfg.flank.mindist);
          end

          % the flanks are shifted by one sample due to the diff
          peaks = peaks + 1;
          % the begin of each waveform is shifted, to ensure that the spike is in the middle
          peaks = peaks - cfg.flank.offset;
          % ensure that no spikes are found within the first and last segment of the data
          % since it is not possible to extract complete waveforms there
          peaks(find(peaks<32))               = [];
          peaks(find((length(dat)-peaks)<32)) = [];

          % --- store the thres val for current channnel
          cfg.chanvals(find(chansel==i),1:3) = [cfg.flank.value NaN cfg.flank.mindist];

      end % cfg.method

      fprintf('detected %d (avg. rate: %.2f)', length(peaks),  (length(peaks) / (length(dat)/hdr.Fs) ));

      % check that n detected spikes is "reasonable" (more than 5
      % percent, less than 80%), otherwise changes the settings
      % note: 32 samples per waveform, length(dat)/32 possible waveforms, hdr.Fs
      if ~strcmp(cfg.adjustselection, 'yes')
        % no automatic adjustment needs to be done
        break;
      elseif numadjustment>10
        % do not auto-adjust more than 10 times
        break;
      else
        adjustValue = [];
        if ( (length(peaks) / (length(dat)/hdr.Fs) )  < 4)
          fprintf(', less than avg. rate of 4 spikes per sec. detected.\n');
          adjustValue = 1+(numadjustment*0.1);
        elseif ~strcmp(cfg.method,'all') && ( (length(peaks) / (length(dat)/hdr.Fs) )  > 600)
          fprintf(', more than avg. rate of 600 spikes per sec. detected.\n');
          adjustValue = 1-(numadjustment*0.1);
        else
          % the detected spike rate is "reasonable", no further adjustments neccessary
          break;
        end

        if ~isempty(adjustValue)
          maxDat = max(abs(zdat))*0.95; % the minimum threshold value to ensure thres. is in correct range
          if strcmp(cfg.method,'zthresh')
            cfg.zthresh.neg = max([cfg.zthresh.neg*adjustValue  -maxDat]);
            cfg.zthresh.pos = min([cfg.zthresh.pos*adjustValue   maxDat]);
            fprintf('... adjusted thresh. to %.2f / %.2f\n',cfg.zthresh.neg, cfg.zthresh.pos);
          elseif strcmp(cfg.method,'flank')
            cfg.flank.value = min([cfg.flank.value*adjustValue  maxDat]);
            fprintf('... adjusted thresh. to %.2f\n',cfg.flank.value);
          end
        end
        numadjustment = numadjustment + 1;
      end

    end % while automatic threshold adjustment

    if strcmp(cfg.interactive, 'no')
      % construct a structure like this
      %   spike.label     = 1xNchans cell-array, with channel labels
      %   spike.waveform  = 1xNchans cell-array, each element contains a matrix (Nsamples X Nspikes), can be empty
      %   spike.timestamp = 1xNchans cell-array, each element contains a vector (1 X Nspikes)
      %   spike.unit      = 1xNchans cell-array, each element contains a vector (1 X Nspikes)
      % note that it only contains a single channel

      spike            = [];
      if isempty(cfg.channelprefix)
        % the label should be a cell-array of length one
        spike.label     = hdr.label(i); 
      else
        % add a prefix to the channel name
        spike.label     = {[cfg.channelprefix '_' hdr.label{i}]};
      end
      spike.waveform   = {zeros(32,length(peaks))};  % FIXME implement variable length waveforms
      spike.timestamp  = {zeros(1,length(peaks), class(hdr.FirstTimeStamp))};
      spike.unit       = {zeros(1,length(peaks))};

      for j=1:length(peaks)
        begsmp = peaks(j);
        endsmp = peaks(j) + 32 - 1;   % FIXME implement a peak shift
        spike.waveform{1}(:,j) = dat(begsmp:endsmp);
        spike.timestamp{1}(j)  = hdr.FirstTimeStamp + typecast((peaks(j)-1)*hdr.TimeStampPerSample, class(hdr.FirstTimeStamp));
      end

      % write the spike data to a new file
      datafile = fullfile(cfg.output, spike.label{1});  % this is without filename extension
      fprintf(', writing to %s\n', datafile);
      ft_write_spike(datafile, spike, 'dataformat', cfg.dataformat, 'fsample', hdr.Fs, 'TimeStampPerSample', hdr.TimeStampPerSample*hdr.Fs);

      % jump out of the interactive loop
      runloop = false;
      newdata = false;

    elseif strcmp(cfg.interactive, 'yes')
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % show the spike data on screen
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      clf
      fprintf('\n\n');
      fprintf('------------------------------------------------------------\n');
      fprintf('Channel %s, detected spike rate is %d per second, ', hdr.label{i}, round(length(peaks)./(cfg.latency(2)-cfg.latency(1))));
      fprintf('change settings or press return to accept\n');
      fprintf('------------------------------------------------------------\n');

      % use the z-transformed data
      dat = zdat;

      subplot('position', [0.1 0.3 0.1 0.6]);
      [hdat, hval] = hist(double(dat), 256);
      h = plot(hdat, hval);
      x = [0 max(hdat)*2];
      if strcmp(cfg.method, 'zthresh')
        y = [cfg.zthresh.neg cfg.zthresh.neg];
        line(x, y, 'color', 'g');
        y = [cfg.zthresh.pos cfg.zthresh.pos];
        line(x, y, 'color', 'g');
      elseif strcmp(cfg.method, 'flank')
        y = [cfg.flank.value cfg.flank.value];
        line(x, y, 'color', 'g');
      end
      abc = axis; abc(1) = 0; abc(2) = 1.2*max(hdat); axis(abc);

      subplot('position', [0.3 0.3 0.6 0.6]);

      if size(cfg.latency,1)==1
        % create a time axis that matches the data
        time = linspace(cfg.latency(1,1), cfg.latency(1,2), length(dat));
      else
        % the data consiss of multiple concatenated segments, create a dummy time axis
        warning('the displayed time axis does not represent real time in the recording');
        time = (0:(length(dat)-1))/hdr.Fs;
      end
      h = plot(time, dat);

      if strcmp(cfg.method, 'zthresh')
        x = [time(1) time(end)];
        y = [cfg.zthresh.neg cfg.zthresh.neg];
        line(x, y, 'color', 'g');
        x = [time(1) time(end)];
        y = [cfg.zthresh.pos cfg.zthresh.pos];
        line(x, y, 'color', 'g');
      elseif strcmp(cfg.method, 'flank')
        x = [time(1) time(end)];
        y = [cfg.flank.value cfg.flank.value];
        line(x, y, 'color', 'g');
      end

      spiketime = (peaks-1)/hdr.Fs;
      spiketime = spiketime + time(1);
      hold on
      plot(spiketime, zeros(size(peaks)), 'r.');
      hold off
      abc = axis; abc(1) = time(1); abc(2) = time(end); axis(abc);

      subplot('position', [0.3 0.1 0.6 0.1]);
      plot(spiketime, zeros(size(peaks)), '.');
      abc = axis; abc(1) = time(1); abc(2) = time(end); axis(abc);

      % start with a wider figure already, this prevents manual resizing
      set(gcf,'Position',[42   302   879   372],'Color','w')

      % ask for a new latency selection
      tmp=(cfg.latency'); oldval = sprintf('%.1f %.1f, ',tmp(:));
      [cfg.latency, newdata] = smartinput(['cfg.latency [' oldval '] = '], cfg.latency);
      runloop = newdata;

      fn = fieldnames(getfield(cfg, cfg.method));
      for k=1:length(fn)
        eval(['oldval = cfg.' cfg.method '.' fn{k} ';']);
        if isnumeric(oldval), oldvalinf = sprintf('%.1f',oldval); else, oldvalinf = oldval; end
        [newval, changed] = smartinput(sprintf('cfg.%s.%s [%s]= ', cfg.method, fn{k},oldvalinf), oldval);
        eval(['cfg.' cfg.method '.' fn{k} ' = newval;']);
        runloop = (runloop | changed);
      end

      if ~runloop
        warning('detected spikes are not written to disk in interactive mode');
      end
      
    end % elseif interactive
    
  end % while runloop

end % for each file

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble provenance

