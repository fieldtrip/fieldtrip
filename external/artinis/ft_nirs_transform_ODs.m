function [data] = ft_nirs_transform_ODs(cfg, data)

% FT_NIRS_TRANSFORM_ODs computes the transformation from optical densities (OD)
% to chromophore concentrations such as (de-)oxygenated  hemoglobin, or the
% other way around.
%
% Use as either
%   [data]      = ft_nirs_transform_ODs(cfg, data)
%   [freq]      = ft_nirs_transform_ODs(cfg, freq)
%   [timelock]  = ft_nirs_transform_ODs(cfg, timelock)
%   [component] = ft_nirs_transform_ODs(cfg, component)
%
%  The configuration "cfg" is a structure containing information about
%  target of the transformation. The configuration should contain
%  cfg.channel            = Nx1 cell-array with selection of channels
%                           (default = 'nirs'), see FT_CHANNELSELECTION for
%                           more details
%  cfg.target             = Mx1 cell-array, can be 'O2Hb' (oxygenated hemo-
%                           globin), 'HHb' de-oxygenated hemoglobin') or
%                           'tHb' (total hemoglobin), or a combination of
%                           those (default: {'O2Hb', 'HHb'})
%
% Optional configuration settings are
%  cfg.age                = scalar, age of the subject (necessary to
%                           automatically select the appropriate DPF, or
%  cfg.dpf                = scalar, differential path length factor
%  cfg.dpffile            = string, location to a lookup table for the
%                           relation between participant age and DPF
%
% Note that the DPF might be different across channels, and is usually
% stored in the optode structure contained in the data.
%
% The function returns data transformed to the specified chromophore
% concentrations of the the requested
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_NIRS_PREPARE_ODTRANSFORMATION, FT_APPLY_MONTAGE

% Options to be implemented:
%   cfg.opto          = structure with optode definition or filename, see FT_READ_SENS
%   cfg.siunits       = yes/no, ensure that SI units are used consistently
%   cfg.logarithm     = string, can be 'natural' or 'base10' (default = 'base10')
%   cfg.inverse       = string, whether the multiplicative inverse should be take
%
% Note that HomER uses the inverse, natural logarithm. The inverse is taken
% after the ln transform.

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
%     Share � copy and redistribute the material in any medium or format
%     Adapt � remix, transform, and build upon the material
%     for any purpose, even commercially.
%
%     The licensor cannot revoke these freedoms as long as you follow the
%     license terms.
%
% Under the following terms:
%
%     Attribution � You must give appropriate credit, provide a link to
%                    the license, and indicate if changes were made. You
%                    may do so in any reasonable manner, but not in any way
%                    that suggests the licensor endorses you or your use.
%
%     ShareAlike � If you remix, transform, or build upon the material,
%                   you must distribute your contributions under the same
%                   license as the original.
%
%     No additional restrictions � You may not apply legal terms or
%                                   technological measures that legally
%                                   restrict others from doing anything the
%                                   license permits.
%
% -----------------------------------
%
% This toolbox is not to be used for medical or clinical purposes.
%
% Copyright (c) 2015-2016 by Artinis Medical Systems.
% Contact: askforinfo@artinis.com
%
% Main programmer: J�rn M. Horschig
% $Id$

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble loadvar     data
ft_preamble provenance  data
ft_preamble trackconfig
ft_preamble debug

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

data = ft_checkdata(data, 'datatype', {'raw', 'timelock', 'freq', 'comp'}, 'feedback', 'yes');

cfg.channel = ft_getopt(cfg, 'channel', 'nirs');

% determine the type of input data
if isfield(data, 'trial')
  israw      = 1;
  isfreq     = 0;
  iscomp     = 0;
  istimelock = 0;
elseif isfield(data, 'freq')
  israw      = 0;
  isfreq     = 1;
  iscomp     = 0;
  istimelock = 0;
elseif isfield(data, 'topo')
  israw      = 0;
  isfreq     = 0;
  iscomp     = 1;
  istimelock = 0;
elseif isfield(data, 'avg')
  israw      = 0;
  iscomp     = 0;
  isfreq     = 0;
  istimelock = 1;
else
  error('input data is not recognized');
end

data = ft_selectdata(cfg, data);

cfg.target  = ft_getopt(cfg, 'target', {'O2Hb', 'HHb'});

target = cfg.target;
if ~iscell(target)
  target     = {target};
end

computetHb  = any(ismember(target, 'tHb'));
computeHHb  = any(ismember(target, 'HHb'));
computeO2Hb = any(ismember(target, 'O2Hb'));

% cfg-handling is done inside here
[montage, cfg] = ft_nirs_prepare_ODtransformation(cfg, data);

% save montage in the cfg
cfg.montage = montage;

% apply the (combined) montages
dataout = ft_apply_montage(data, montage);

if computetHb
  % total hemoglobin is the sum of oxygenated and deoxygenated hemoglobin
  tmpcfg = [];
  tmpcfg.channel = '*[O2Hb]';
  dataO2Hb = ft_selectdata(tmpcfg, dataout);
  dataO2Hb.label = strrep(dataO2Hb.label, 'O2Hb', 'tHb');

  tmpcfg = [];
  tmpcfg.channel = '*[HHb]';
  dataHHb = ft_selectdata(tmpcfg, dataout);
  dataHHb.label = strrep(dataHHb.label, 'HHb', 'tHb');

  tmpcfg = [];
  tmpcfg.operation = 'add';
  tmpcfg.parameter = 'trial';
  dataTotal = ft_math(tmpcfg, dataO2Hb, dataHHb);

  if ~computeHHb && ~computeO2Hb
    % replace dataout
    dataout.label = dataTotal.label;
    dataout.trial = dataTotal.trial;
  else
    % concat to dataout
    dataout.label = vertcat(dataout.label, dataTotal.label);
    dataout.trial = cellfun(@vertcat, dataout.trial, dataTotal.trial, 'UniformOutput', false);
  end
end

% remove O2 or H if not desired
tmpcfg = [];
tmpcfg.channel = [];
if computetHb
  tmpcfg.channel = horzcat(tmpcfg.channel, {'*tHb]'});
end
if computeHHb
  tmpcfg.channel =  horzcat(tmpcfg.channel, {'*HHb]'});
end
if computeO2Hb
  tmpcfg.channel =  horzcat(tmpcfg.channel, {'*O2Hb]'});
end

dataout = ft_selectdata(tmpcfg, dataout);

if size(dataout)>1
  if isfreq
    data = ft_appendfreq([], dataout{:});
  elseif istimelock
    data = ft_appendtimelock([], dataout{:});
  else
    data = ft_appenddata([], dataout{:});
  end
else
  data = dataout;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous    data
ft_postamble provenance  data
ft_postamble history     data
