function [data] = ft_transform_ODs(cfg, data)

% FT_TRANSFORM_ODs computes the transformation from
% optical densities (OD) to chromophore concentrations such as (de-)
% oxygenated hemoglobin, or the other way around.
%
% Use as either
%   [data] = ft_transform_ODs(cfg, data);
%   [freq] = ft_transform_ODs(cfg, freq);
%   [timelock] = ft_transform_ODs(cfg, timelock);
%   [component] = ft_transform_ODs(cfg, component);
%
%  The configuration "cfg" is a structure containing information about
%  target of the transformation. The configuration should contain
%  cfg.channel            = Nx1 cell-array with selection of channels
%                           (default = 'nirs'), see FT_CHANNELSELECTION for
%                           more details
%  cfg.target             = Mx1 cell-array, can be 'OD' (optical densities), 
%                           'O2Hb' (oxygenated hemoglobin), 'HHb' de-
%                           oxygenated hemoglobin') or 'tHb' (total hemo-
%                           globin), or a combination of those 
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
% See also FT_PREPARE_ODTRANSFORMATION, FT_APPLY_MONTAGE

% options to be implemented:
% 
% The NIRS positions can be present in the data or can be specified as
%   cfg.opto          = structure with optode positions, see FT_DATATYPE_SENS
%   cfg.optofile      = name of file containing the optode positions, see FT_READ_SENS
%
%   cfg.siunits       = ft_getopt(cfg, 'siunits', 'no');   % yes/no, ensure that SI units are used consistently
%
%   cfg.logarithm     = string, can be 'natural' or 'base10' (default = 'base10')
%   cfg.inverse       = string, whether the multiplicative inverse should be take
%
% Note that HomER uses the inverse, natural logarithm. The inverse is taken after the ln transform

% This file is part of the FieldTrip NIRS toolbox developed and maintained
% by Artinis Medical Systems (http://www.artinis.com). For more information
% on FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-NoDerivs 3.0 Unported License. To view a copy
% of this license, visit http://creativecommons.org/licenses/by-nc-nd/3.0/
% or send a letter to Creative Commons, PO Box 1866, Mountain View, CA
% 94042, USA.
% 
% Creative Commons Attribution-NonCommercial-NoDerivs 3.0 Unported License:
% -----------------------------------
% You are free to:
%     Share — copy and redistribute the material in any medium or format
%     The licensor cannot revoke these freedoms as long as you follow the license terms.
% 
% Under the following terms:
%     Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
%     NonCommercial — You may not use the material for commercial purposes.
%     NoDerivatives — If you remix, transform, or build upon the material, you may not distribute the modified material.
%     No additional restrictions — You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.
% -----------------------------------
% 
% This toolbox is not to be used for medical or clinical purposes.
% 
% Copyright (c) 2015 by Artinis Medical Systems.
% Contact: askforinfo@artinis.com
%
% Main programmer: Jörn M. Horschig
% $Id: ft_transform_ODs.m 10169 2015-02-05 12:56:47Z jorhor $

revision = '$Id: ft_transform_ODs.m 10169 2015-02-05 12:56:47Z jorhor $';


% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble loadvar     data
ft_preamble provenance  data
ft_preamble trackconfig
ft_preamble debug

data = ft_checkdata(data, 'datatype', {'raw', 'timelock', 'freq', 'comp'}, 'feedback', 'yes');

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

cfg.target  = ft_getopt(cfg, 'target', 'tHb');

target = cfg.target;
if ~iscell(target)
  target     = {target};
end

tHbidx = ismember(target, 'tHb');
HHbidx = ismember(target, 'HHb');
O2Hbidx = ismember(target, 'O2Hb');

computetHb = tHbidx~=0;
existsHHb  = isempty(HHbidx);
existsO2Hb  = isempty(O2Hbidx);

% total hemoglobin is the sum of oxygenated and deoxygenated hemoglobin
if computetHb
  if ~existsHHb
    target{end+1} = 'HHb';
    HHbidx = numel(target);
  end
  if ~existsO2Hb
    target{end+1} = 'O2Hb';
    O2Hbidx = numel(target);
  end
end


for t=1:numel(target)
  
  if (t==tHbidx)
    continue; % this will be dealt with later
  end
  
  cfg.target = target{t};
  % cfg-handling is done inside here
  [montage, cfg] = ft_prepare_ODtransformation(cfg, data);
  
  % save montage in the cfg
  cfg.montage{t} = montage;  
end

if computetHb
  % make a single montage matrix
  cfg.montage{tHbidx}.labelorg = cfg.montage{HHbidx}.labelorg;
  cfg.montage{tHbidx}.labelnew = regexprep(cfg.montage{HHbidx}.labelnew, sprintf('%s%s', regexptranslate('wildcard', '[*]'), '$'), sprintf('[%s]', target{tHbidx}));
  cfg.montage{tHbidx}.tra      = cfg.montage{HHbidx}.tra+cfg.montage{O2Hbidx}.tra;  
  
  %order reversed wrt above!
  if ~existsO2Hb
    cfg.montage(O2Hbidx) = [];
    target(O2Hbidx) = [];
  end
  
  if ~existsHHb
    cfg.montage(HHbidx) = [];
    target(HHbidx) = [];
  end
end

dataout = [];
cfg.target = target;
for t=1:numel(cfg.target)
  % apply the (combined) montages
  dataout{t} = ft_apply_montage(data, cfg.montage{t});
end

if size(dataout)>1
  if isfreq
    data = ft_appendfreq([], dataout{:});
  elseif istimelock
    data = ft_appendtimelock([], dataout{:});
  else
    data = ft_appenddata([], dataout{:});
  end
else
  data = dataout{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous    data
ft_postamble provenance  data
ft_postamble history     data