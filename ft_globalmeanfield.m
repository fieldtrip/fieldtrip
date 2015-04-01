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
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

revision = '$Id$';

% do the general setup of the function

% the ft_preamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_defaults                 % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init            % this will reset warning_once and show the function help if nargin==0 and return an error
ft_preamble provenance      % this records the time and memory usage at the beginning of the function
ft_preamble trackconfig     % this converts the cfg structure in a config object, which tracks the cfg options that are being used
ft_preamble debug           % this allows for displaying or saving the function name and input arguments upon an error
ft_preamble loadvar datain  % this reads the input data in case the user specified the cfg.inputfile option

% the abort variable is set to true or false in ft_preamble_init
if abort
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
method    = ft_getopt(cfg, 'method', 'amplitude');        

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

% the ft_postamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_postamble debug            % this clears the onCleanup function used for debugging in case of an error
ft_postamble trackconfig      % this converts the config object back into a struct and can report on the unused fields
ft_postamble provenance       % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and MATLAB version etc. to the output cfg
ft_postamble previous datain  % this copies the datain.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble history dataout  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
ft_postamble savevar dataout  % this saves the output data structure to disk in case the user specified the cfg.outputfile option
