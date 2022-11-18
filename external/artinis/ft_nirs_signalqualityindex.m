function dataout = ft_nirs_signalqualityindex(cfg, datain)

% FT_NIRS_SIGNALQUALITYINDEX processes NIRS data and returns a copy of the original
% data replaced with nans in the signal segements that are bellow the specified
% quality threshold.
%
% Use as
%   dataout = ft_nirs_signalqualityindex(cfg, datain)
% where cfg is a configuration structure and indata is raw NIRS-data (in optical
% densities, ODs) that is represented according to the output of FT_PREPROCESSING.
%
% The configuration should contain the following options
%   cfg.threshold    = scalar, the SQI (signal quality index) value that
%                      has to be exceeded to be labelled as a 'good'
%                      channel (default = 3.5)
%   cfg.windowlength = scalar, the length (in seconds) of the signal
%                      segments to be analyzed (default = 10)
%
% This function is based on:
% - Sappia, M. S., Hakimi, N., Colier, W. N., & Horschig, J. M. (2020).
%   Signal quality index: an algorithm for quantitative assessment of
%   functional near infrared spectroscopy signal quality. Biomedical Optics
%   Express, 11(11), 6732-6754. https://doi.org/10.1364/BOE.409317
%
% Please cite accordingly. Thank you!
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the% input/output structure.
%
% See also FT_NIRS_SCALPCOUPLINGINDEX, FT_NIRS_REFERENCECHANNELSUBTRACTION, FT_NIRS_TRANSFORM_ODS

% You are using the FieldTrip NIRS toolbox developed and maintained by
% Artinis Medical Systems (http://www.artinis.com). For more information
% on FieldTrip, see http://www.fieldtriptoolbox.org
%
% This work is licensed under a Creative Commons Attribution-ShareAlike 4.0
% International License. To view a copy of this license, visit
% http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
%
% Creative Commons Attribution-ShareAlike 4.0 International License:
% -----------------------------------
% You are free to:
%
%     Share - copy and redistribute the material in any medium or format
%     Adapt - remix, transform, and build upon the material
%     for any purpose, even commercially.
%
%     The licensor cannot revoke these freedoms as long as you follow the
%     license terms.
%
% Under the following terms:
%
%     Attribution - You must give appropriate credit, provide a link to
%                    the license, and indicate if changes were made. You
%                    may do so in any reasonable manner, but not in any way
%                    that suggests the licensor endorses you or your use.
%
%     ShareAlike - If you remix, transform, or build upon the material,
%                   you must distribute your contributions under the same
%                   license as the original.
%
%     No additional restrictions - You may not apply legal terms or
%                                   technological measures that legally
%                                   restrict others from doing anything the
%                                   license permits.
%
% -----------------------------------
%
% This toolbox is not to be used for medical or clinical purposes.
%
% Copyright (c) 2021 by Artinis Medical Systems.
% Contact: askforinfo@artinis.com
%
% Main programmers:
% Naser Hakimi, Artinis Medical Systems BV, http://www.artinis.com
% Sofia Sappia, Artinis Medical Systems BV, http://www.artinis.com
% $Id$


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function

% the ft_preamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    datain
ft_preamble provenance datain


% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% ensure that the input data is raw NIRS-data, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
datain = ft_checkdata(datain, 'datatype', 'raw', 'senstype', 'nirs');

% get the options
cfg.threshold    = ft_getopt(cfg, 'threshold', 3.5);
cfg.windowlength = ft_getopt(cfg, 'windowlength', 10);
% ensure that the options are valid
cfg = ft_checkopt(cfg, 'threshold', 'doublescalar');
cfg = ft_checkopt(cfg, 'windowlength', 'doublescalar');

data_resampled = datain; % copy of data to compute SQI on
dataout = datain; % where low quality signal segments will be replaced with nans

% check if labels match the expected order (wavelengths belonging to the
% same channel are consecutive)
matching_channels = match_nirs_labels(data_resampled.label);
idx_channels = 1:size(matching_channels,1);

for idx_temp_channel=1:size(matching_channels,1)
  idx_matching = find(matching_channels(idx_temp_channel,:)==1);
  if numel(idx_matching) ~= 2
    error('Expected channel combination associated with 2 wavelengths. Found %d',numel(idx_matching));
  end
  idx_channel1 = idx_matching(1);
  idx_channel2 = idx_matching(2);
  if abs(idx_channel2 - idx_channel1) > 1
    error('Labels do not match the expected order');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                signal quality index computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if sample rate is not 50 Hz, then resample
if data_resampled.fsample ~= 50
  res_cfg = [];
  res_cfg.resamplefs = 50;
  data_resampled = ft_resampledata(res_cfg, data_resampled);
end

% convert to oxy and deoxy hemoglobin
cfg_Hb = [];
[data_Hb] = ft_nirs_transform_ODs(cfg_Hb, data_resampled);

resampled_windowsize = cfg.windowlength * data_resampled.fsample; % in samples
org_windowsize = cfg.windowlength * datain.fsample;

% loop through all channels and trials
num_channels = numel(data_resampled.label);
num_trials = numel(data_resampled.trial);

while(numel(idx_channels)>0)
  idx_temp_channel=idx_channels(1);
  idx_matching = find(matching_channels(idx_temp_channel,:)==1);
  idx_channel1 = idx_matching(1);
  idx_channel2 = idx_matching(2);
  idx_channels(idx_channels == idx_matching(1)) = [];
  idx_channels(idx_channels == idx_matching(2)) = [];
  for idx_trial = 1:num_trials
    % get ODs and Hb signals per trial, per channel
    OD1 = data_resampled.trial{idx_trial}(idx_channel1,:);
    OD2 = data_resampled.trial{idx_trial}(idx_channel2,:);
    oxy = data_Hb.trial{idx_trial}(idx_channel1,:);
    dxy = data_Hb.trial{idx_trial}(idx_channel2,:);
    
    % compute window onset samples per signal segment
    window_onsets_resampled = 1:resampled_windowsize:numel(OD1) - resampled_windowsize; % no overlap
    window_onsets_org = 1:org_windowsize:size(dataout.trial{idx_trial}(idx_channel1,:),2) - org_windowsize;
    
    for idx_window = 1:numel(window_onsets_resampled)
      idx_signal_segment = window_onsets_resampled(idx_window):window_onsets_resampled(idx_window) + resampled_windowsize;
      % get signal segments
      OD1_segment = OD1(idx_signal_segment);
      OD2_segment = OD2(idx_signal_segment);
      oxy_segment = oxy(idx_signal_segment);
      dxy_segment = dxy(idx_signal_segment);
      % compute SQI score
      sqi_score = SQI(OD1_segment, OD2_segment, oxy_segment, dxy_segment, data_resampled.fsample);
      if sqi_score < cfg.threshold
        idx_sample_start = window_onsets_org(idx_window);
        idx_signal_segment_orig = idx_sample_start:idx_sample_start + org_windowsize;
        dataout.trial{idx_trial}(idx_channel1, idx_signal_segment_orig) = NaN;
        dataout.trial{idx_trial}(idx_channel2, idx_signal_segment_orig) = NaN;
      end
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            end of signal quality index computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the ft_postamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_postamble debug

ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout
