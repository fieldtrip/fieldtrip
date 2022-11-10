function dataout = ft_nirs_scalpcouplingindex(cfg, datain)

% FT_NIRS_SCALPCOUPLINGINDEX computes the zero-lag cross-correlation between pairs of
% raw NIRS-channels to identify bad channels. In the output data structure the bad
% channels are removed, or filled with NaNs.
%
% Use as
%   [outdata] = ft_scalpcouplingindex(cfg, indata)
% where cfg is a configuration structure and indata is raw NIRS-data (in optical
% densities, ODs) that is represented according to the output of FT_PREPROCESSING.
%
% The configuration should contain the following options
%   cfg.threshold   = scalar, the correlation value which has to be exceeded to be
%                     labelled a 'good' channel (default = 0.75)
%   cfg.keepchannel = string, determines how to deal with channels that are not selected, can be
%                     'no'  completely remove deselected channels from the data (default)
%                     'nan' fill the channels that are deselected with NaNs
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% This function is based on:
% - Pollonini, L., Olds, C., Abaya, H., Bortfeld, H., Beauchamp, M. S., & Oghalai, J. S. (2014).
%   Auditory cortex activation to natural speech and simulated cochlear implant speech measured
%   with functional near-infrared spectroscopy. Hearing Research, 309, 84?93. doi:10.1016/j.heares.2013.11.007
%
% Please cite accordingly. Thank you!
%
% See also FT_NIRS_SIGNALQUALITYINDEX, FT_NIRS_REFERENCECHANNELSUBTRACTION, FT_NIRS_TRANSFORM_ODS

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
% Copyright (c) 2016 by Artinis Medical Systems.
% Contact: askforinfo@artinis.com
%
% Main programmer:
% Jörn M. Horschig, Artinis Medical Systems BV, http://www.artinis.com
% Marc van Wanrooij, DCN, http://www.neural-code.com
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function

% the ft_preamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_defaults                   % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init              % this will reset ft_warning and show the function help if nargin==0 and return an error
ft_preamble debug             % this allows for displaying or saving the function name and input arguments upon an error
ft_preamble loadvar    datain % this reads the input data in case the user specified the cfg.inputfile option
ft_preamble provenance datain % this records the time and memory usage at the beginning of the function

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
cfg.threshold    = ft_getopt(cfg, 'threshold', 0.75);
cfg.keepchannel = ft_getopt(cfg, 'keepchannel', 'no');

% ensure that the options are valid
cfg = ft_checkopt(cfg, 'threshold', 'doublescalar');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if sample rate is too high (>10Hz) or we hit Nyquist then resample
dataf = datain;
if dataf.fsample > 10 || dataf.fsample <= 5
  rescfg = [];
  rescfg.resamplefs = 10;
  dataf = ft_resampledata(rescfg, dataf);
end

% filter the signal in the range of interest
filcfg = [];
filcfg.bpfilter = 'yes';
filcfg.bpfreq = [0.5 2.5];
dataf = ft_preprocessing(filcfg, dataf);

% loop through all channels, find similar ones and correlate
nTrials = numel(dataf.trial);
nChans = numel(dataf.label);
i = 1;
skipChan = false(nChans);
sci = ones(nChans, 1);
while i<nChans
  % if this is a non-raw-nirs channel, then skip
  if ~strcmp(dataf.label{i}(end-2:end), 'nm]') % all raw nirs channel end with nm]
    skipChan(i) = true;
  end
  
  if skipChan(i)
    i = i+1;
  else
    % find all matching channels
    bracket_offset = strfind(dataf.label{i}, '[');
    chanidx = find(~cellfun(@isempty, strfind(dataf.label, dataf.label{i}(1:bracket_offset-1))));
    % correlate for each pair
    for ca=1:numel(chanidx)-1
      cidxA = chanidx(ca);
      if skipChan(cidxA)
        continue;
      end
      
      % initialize this channel
      sci(cidxA) = 0;
      for cb = ca+1:numel(chanidx)
        cidxB = chanidx(cb);
        if skipChan(cidxB)
          continue;
        end
        
        % initialize this channel
        sci(cidxB) = 0;
        
        % correlate for each trial
        rsum = 0;
        for tr=1:nTrials
          r		= corrcoef(dataf.trial{tr}(cidxA, :), dataf.trial{tr}(cidxB, :));
          rsum		= rsum + r(2);
        end
        rsum = rsum/nTrials;
        
        % sum up the sci, divided by number of 'companion wavelengths'
        sci(cidxA) = sci(cidxA) + rsum/(numel(chanidx)-1);
        sci(cidxB) = sci(cidxB) + rsum/(numel(chanidx)-1);
      end % end for:cb
      
      % mark this channel as 'done'
      skipChan(cidxA) = true;
    end % end for:ca
    
    % skip the last channel as it's either 'done' or has no 'companion'
    skipChan(chanidx(end)) = true;
    
  end % end if:skipChan
  
end

% remove the bad channels on the original datain (default) or replaces with
% nans
chanidx = sci > cfg.threshold;
goodchannel = datain.label( chanidx);
badchannel  = datain.label(~chanidx);

if ~isempty(badchannel)
  switch cfg.keepchannel
    case 'no'
      selcfg = [];
      selcfg.channel = goodchannel;
      dataout = ft_selectdata(selcfg, datain);
      fprintf('the following channels were removed: '); % to be continued below ...
    case 'nan'
      dataout = datain;
      for i = 1:length(datain.trial)
        dataout.trial{i}(~chanidx,:) = nan;
      end
      fprintf('the following channels were filled with NaNs: '); % to be continued below ...
    otherwise
      error('invalid option for cfg.keepchannel');
  end
  
  % show which channels were removed or filled with NaNs
  for i=1:(length(badchannel)-1)
    fprintf('\n%s', badchannel{i});
  end
  fprintf('\n%s\n', badchannel{end});
  
else
  dataout=datain;
end

% this might involve more active checking of whether the input options
% are consistent with the data and with each other

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function

% the ft_postamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_postamble debug               % this clears the onCleanup function used for debugging in case of an error
ft_postamble previous   datain   % this copies the datain.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble provenance dataout  % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and MATLAB version etc. to the output cfg
ft_postamble history    dataout  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
ft_postamble savevar    dataout  % this saves the output data structure to disk in case the user specified the cfg.outputfile option
