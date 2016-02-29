function dataout = ft_globalmeanfield(cfg, datain)

% FT_GLOBALMEANFIELD calculates global mean field amplitude or power of input data
%
% Use as
%   [gmf] = ft_globalmeanfield(cfg, indata)
%
% The data should be organised in a structure as obtained from the
% FT_TIMELOCKANALYSIS function. The configuration should be according to
% FT_PREPROCESSING function. The configuration should be according to
%
%   cfg.methods   = string, whether the amplitude or power should be calculated (see below, default is 'amplitude')
%   cfg.channel   = Nx1 cell-array with selection of channels (default = 'all'),
%                            see FT_CHANNELSELECTION for details
%
% This function calculates the global mean field power, or amplitude,
% as described in:
% Lehmann D, Skrandies W. Reference-free identification of components of
% checkerboard-evoked multichannel potential fields. Electroencephalogr Clin
% Neurophysiol. 1980 Jun;48(6):609-21. PubMed PMID: 6155251.
%
% Please note that to calculate what is clasically referred to as Global
% Mean Field Power, cfg.method must be 'amplitude'. The naming implies a
% squared measure but this is not the case.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_TIMELOCKANALYSIS

% Copyright (C) 2014, Jim Herring
% Copyright (C) 2014, Robert Oostenveld
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar datain
ft_preamble provenance datain
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'channel'});
datain = ft_selectdata(tmpcfg, datain);
[cfg, datain] = rollback_provenance(cfg, datain);

% ensure that the input data is valid for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
datain = ft_checkdata(datain, 'datatype', {'timelock'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% get the options
cfg.method    = ft_getopt(cfg, 'method', 'amplitude');

% ensure that the options are valid
cfg = ft_checkopt(cfg, 'method', 'char', {'amplitude', 'power'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% do your stuff...
dataout = [];
dataout = keepfields(datain,{'avg','time','label','dimord','cfg'});

switch cfg.method
  case 'amplitude'
    dataout.avg = std(dataout.avg,1);
  case 'power'
    dataout.avg = std(dataout.avg,1).^2;
end;

dataout.label = {'gmfp'};

% this might involve more active checking of whether the input options
% are consistent with the data and with each other

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout
