function dataout = ft_examplefunction(cfg, datain)

% FT_EXAMPLEFUNCTION demonstrates to new developers how a FieldTrip function should look like
%
% Use as
%  outdata = ft_examplefunction(cfg, indata) 
% where indata is <<describe the type of data or where it comes from>> 
% and cfg is a configuratioun structure that should contain 
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
% Seee also <<give a list of function names, all in capitals>>

% Here come the Copyrights
%
% Here comes the Revision tag, which is auto-updated by the version control system
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_defaults

hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    datain = loadvar(cfg.inputfile, 'data');
  end
end

% ensure that the input data is valiud for this function, this will also do 
% backward-compatibility conversions of old data that for example was 
% read from an old *.mat file
datain = ft_checkdata(datain, 'datatype', {'raw', 'comp'}, 'feedback', 'yes', 'hastrialdef', 'yes', 'hasoffset', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
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

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add the version details of this function call to the configuration
cfg.version.name = mfilename('fullpath'); % this is helpful for debugging
cfg.version.id   = '$Id$'; % this will be auto-updated by the revision control system

% add information about the Matlab version used to the configuration
cfg.version.matlab = version(); % this is helpful for debugging

if hasdata && isfield(data, 'cfg')
  % remember the configuration details of the input data
  cfg.previous = datain.cfg;
end

% remember the exact configuration details in the output
dataout.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', dataout); % use the variable name "data" in the output file
end
