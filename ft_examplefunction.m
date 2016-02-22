function dataout = ft_examplefunction(cfg, datain)

% FT_EXAMPLEFUNCTION demonstrates to new developers how a FieldTrip function should look like
%
% Use as
%   outdata = ft_examplefunction(cfg, indata)
% where indata is <<describe the type of data or where it comes from>>
% and cfg is a configuration structure that should contain
%
% <<note that the cfg list should be indented with two spaces
%
%  cfg.option1    = value, explain the value here (default = something)
%  cfg.option2    = value, describe the value here and if needed
%                   continue here to allow automatic parsing of the help
%
% The configuration can optionally contain
%   cfg.option3   = value, explain it here (default is automatic)
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

% ensure that the input data is valid for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
datain = ft_checkdata(datain, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'deprecated',  {'normalizecov', 'normalizevar'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blcwindow', 'baselinewindow'});

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'method', 'foi', 'tapsmofrq'});

% ensure that the options are valid
cfg = ft_checkopt(cfg, 'vartrllen', 'double', {0, 1, 2});
cfg = ft_checkopt(cfg, 'method', 'char', {'mtm', 'convol'});

% get the options
method    = ft_getopt(cfg, 'method');        % there is no default
vartrllen = ft_getopt(cfg, 'vartrllen', 2);  % the default is 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do your stuff...
dataout = [];

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
