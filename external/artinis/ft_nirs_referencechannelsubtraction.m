function data_rcs = ft_nirs_referencechannelsubtraction(cfg, datain)

% FT_NIRS_REFERENCECHANNELSUBTRACTION performs reference channel subtraction for NIRS data
%
% Use as
%   outdata = ft_nirs_referencechannelsubtraction(cfg, indata)
% where indata is nirs data and cfg is a configuration structure that should contain
%
%   cfg.shortdistance = scalar, below which distance a channel is regarded
%                       as short in cm (default = 1.5)
%   cfg.closedistance = scalar, defines the maximal distance between a
%                       long and a short channel in cm (default = 15).
%                       NOT APPLIED CURRENTLY!
%   cfg.method        = string, 'regstat2', 'QR' or 'OLS' (default = 'QR')
%   cfg.verbose       = boolean, whether text output is desired (default = false)
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_NIRS_SCALPCOUPLINGINDEX, FT_NIRS_TRANSFORM_ODS

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
% Marc van Wanrooij, DCN, http://www.neural-code.com
% Jörn M. Horschig, Artinis Medical Systems BV, http://www.artinis.com
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
ft_preamble trackconfig       % this converts the cfg structure in a config object, which tracks the cfg options that are being used

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
cfg.shortdistance = ft_getopt(cfg, 'shortdistance', 1.5);
cfg.closedistance = ft_getopt(cfg, 'closedistance', 15);
cfg.method        = ft_getopt(cfg, 'method', 'QR');
cfg.verbose       = ft_getopt(cfg, 'verbose', false);

% check the options
cfg = ft_checkopt(cfg, 'method', 'char', {'regstat2', 'QR', 'OLS'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute distances
distances = ft_nirs_optodedistance(datain);
% compute correlation
selshallow	= distances<cfg.shortdistance;
shallowIdx = find(selshallow);
seldeep		= distances>=cfg.shortdistance;
deepIdx = find(seldeep);

data_rcs			 = datain;
data_rcs.label = datain.label(seldeep);
for tr=1:numel(datain.trial)
  dat = datain.trial{tr};
  shallow		= dat(selshallow,:);
  shallow		= bsxfun(@minus,shallow,mean(shallow,2)); % mean detrend
  shallowlabel = datain.label(selshallow);
  
  deep		= dat(seldeep,:);
  deep		= bsxfun(@minus,deep,mean(deep,2)); % mean detrend
  deeplabel = datain.label(seldeep);
  
  time		= datain.time{tr};
  
  %% Reference channel subtraction
  ndeep		= size(deep,1);
  signal		= NaN(size(deep));
  x			= shallow';
  for dpIdx	= 1:ndeep
    y				= deep(dpIdx,:)';
    switch (cfg.method)
      case 'regstat2'
        b				= regstats2(y,x,'linear',{'beta','r'});
        beta    = b.beta;
        res			= b.r;

      case 'QR'
        [Q,R] = qr(x,0);
        beta = R\(Q'*y);
        yhat = x*beta;
        res = y - yhat;
        
      case 'OLS'
        x2 = [repmat(1, size(x, 1), 1) x];
        beta = x2\y;
        yhat = x2*beta;
        res  = y - yhat;
        beta(1) = [];
        
      otherwise % it should never come here as we use ft_checkopt
        error('unrecognized method');
    end
    
    signal(dpIdx,:) = res';

    % sanity check of results
    if cfg.verbose
      fprintf('Found the following meaningful shallow channels for deep channel %s:"', deeplabel{dpIdx});
      shIdx = find(beta>0.5);
      deepPos = datain.opto.chanpos(deepIdx(dpIdx), :);
      shallowPos = datain.opto.chanpos(shallowIdx(shIdx), :);
      posDist = sqrt(sum((shallowPos - repmat(deepPos, size(shallowPos, 1), 1)).^2, 2));
      closeChan = posDist < cfg.closedistance;
      for s=1:numel(shIdx)
        if (closeChan(s))
          fprintf('\n\t%s (d=%.2fcm)', shallowlabel{shIdx(s)}, posDist(s))
        else
          fprintf('\n\t#%s (d=%.2fcm)', shallowlabel{shIdx(s)}, posDist(s))
        end
      end
      fprintf('\n\n')
    end
  end
  %% overwrite
  data_rcs.time{tr}	= time;
  data_rcs.trial{tr}	= signal;
end

%% Auto-correlation
% rshallow	= corrcoef(shallow'); % correlation matrix shallow channels
% rdeep		= corrcoef(deep'); % correlation matrix deep channels
% rsignal		= corrcoef(signal'); % correlation matrix signal channels


% this might involve more active checking of whether the input options
% are consistent with the data and with each other

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function

% the ft_postamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_postamble debug               % this clears the onCleanup function used for debugging in case of an error
ft_postamble trackconfig         % this converts the config object back into a struct and can report on the unused fields
ft_postamble previous   datain   % this copies the datain.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble provenance dataout  % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and MATLAB version etc. to the output cfg
ft_postamble history    dataout  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
ft_postamble savevar    dataout  % this saves the output data structure to disk in case the user specified the cfg.outputfile option
