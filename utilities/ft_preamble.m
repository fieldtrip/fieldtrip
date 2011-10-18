function ft_preamble(caller_cmd, varargin)

% FT_PREAMBLE is a helper function that is included in many of the FieldTrip
% functions and which takes care of some general settings and operations at the
% begin of the function

% ideally this would be a script, because the local variables would then be
% shared with the calling function. Instead, the variables are explicitely
% copied back and forth between the callers workspace and this one.

%% copy all variables from the calling function to the local workspace
caller_var = evalin('caller', 'whos');
for caller_indx=1:length(caller_var)
  eval(sprintf('%s = evalin(''caller'', ''%s'');', caller_var(caller_indx).name, caller_var(caller_indx).name));
end
clear caller_var caller_indx

% these are also useful to have
nargin  = evalin('caller', 'nargin');
nargout = evalin('caller', 'nargout');

%% execute the desired preamble section
switch caller_cmd
  
  case 'defaults'
    % set up the path and the global ft_defaults
    ft_defaults
    
  case 'callinfo'
    % record the start time and memory
    ftFuncTimer = tic();
    ftFuncClock = clock();
    ftFuncMem   = memtic();
    
  case 'trackconfig'
    % most fieldtrip functions should allow for configuration tracking, except for
    % the functions that take a cfg as input and return a cfg as output
    cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
    
  case 'loadvar'
    % load optional given inputfile as data
    hasdata = (nargin>1);
    if isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
      % the input data should be read from file
      if hasdata
        error('cfg.inputfile should not be used in conjunction with giving input data to this function');
      else
        assign(varargin{1}, loadvar(cfg.inputfile, varargin{1}));
      end
    end
    
  otherwise
    error('unknown preamble command "%s"', caller_cmd);
end

% these should not be copied back to the caller
clear caller_cmd nargin nargout varargin varargout

%% copy all variables from the local workspace back to the calling function
this_var = whos;
for this_indx=1:length(this_var)
  assignin('caller', this_var(this_indx).name, eval(this_var(this_indx).name));
end
clear this_var this_indx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION this compares to assignin('this', ...)
% See also http://woodshole.er.usgs.gov/staffpages/cdenham/public_html/snackbar/assign.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function assign(var, val)
assignin('caller', var, val);
