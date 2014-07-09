function dataout = ft_globalmeanfieldpower(cfg, datain)

% FT_GLOBALMEANFIELDPOWER calculates global mean field power of input data
%
% Use as
%   outdata = ft_globalmeanfieldpower(cfg, indata) 
% where indata is timelock averaged data 
% and cfg is a configuration structure that can contain 
%
% The configuration can optionally contain
%   cfg.methods   = string, whether the output should be squared (default is 'amplitude')
%   cfg.channel   = Nx1 cell-array with selection of channels (default = 'all'),
%                            see FT_CHANNELSELECTION for details
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also <<give a list of function names, all in capitals>>

% Here come the Copyrights
%
% Here comes the Revision tag, which is auto-updated by the version control system
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
ft_preamble provenance      % this records the time and memory usage at teh beginning of the function
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
ft_postamble provenance       % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and matlab version etc. to the output cfg
ft_postamble previous datain  % this copies the datain.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble history dataout  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
ft_postamble savevar dataout  % this saves the output data structure to disk in case the user specified the cfg.outputfile option
