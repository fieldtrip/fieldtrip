function ft_postamble(caller_cmd, varargin)

% FT_POSTAMBLE is a helper function that is included in many of the FieldTrip
% functions and which takes care of some general settings and operations at the
% end of the function

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
nargin    = evalin('caller', 'nargin');
nargout   = evalin('caller', 'nargout');
revision  = evalin('caller', 'revision');
stack     = dbstack('-completenames');
stack     = stack(2);

%% execute the desired preamble section
switch caller_cmd
  
  case 'previous'
    % remember the cfg history of the input data structures
    for i=1:length(varargin)
      if isfield(eval(varargin{1}), 'cfg')
        cfg.previous{i} = getfield(eval(varargin{1}), 'cfg');
      else
        cfg.previous{i} = [];
      end
    end
    clear i
    
  case 'trackconfig'
    % accessing this field here is needed for the configuration tracking
    % by accessing it once, it will not be removed from the output cfg
    evalin('caller', 'try, cfg.outputfile; end');
    % get the output cfg
    cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');
    
  case 'callinfo'
    % add information about the Matlab version used to the configuration
    cfg.callinfo.matlab = version();
    
    % add information about the function call to the configuration
    cfg.callinfo.proctime = toc(ftFuncTimer);
    cfg.callinfo.procmem  = memtoc(ftFuncMem);
    cfg.callinfo.calltime = ftFuncClock;
    cfg.callinfo.user     = getusername();
    % add version information to the configuration
    cfg.version.name = stack.file;
    cfg.version.id   = revision;
    
    % give some feedback
    fprintf('the call to "%s" took %d seconds and an estimated %d MB\n', stack.name, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));
    
  case 'savevar'
    % the output data should be saved to a MATLAB file
    if ~isempty(cfg.outputfile)
      savevar(cfg.outputfile, varargin{1}, eval(varargin{1}));
    end
    
  otherwise
    error('unknown postamble command "%s"', caller_cmd);
end

% these should not be copied back to the caller
clear caller_cmd nargin nargout varargin varargout revision stack

%% copy all variables from the local workspace back to the calling function
this_var = whos;
for this_indx=1:length(this_var)
  assignin('caller', this_var(this_indx).name, eval(this_var(this_indx).name));
end
clear this_var this_indx
