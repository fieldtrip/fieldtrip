% FT_POSTAMBLE_CALLINFO is a helper script that reports the time and memory
% that was used by the calling function and that stores this information
% together with user name, MATLAB version and some other stuff in the
% output cfg structure. This should be used together with FT_PREAMBLE_CALLINFO,
% which records the time and memory at the start of the function.

stack = dbstack('-completenames');
% stack(1) is this script
% stack(2) is the calling ft_postamble function
% stack(3) is the main FieldTrip function that we are interested in
stack = stack(3);

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();

% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.procmem  = memtoc(ftFuncMem);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user     = getusername();
cfg.callinfo.hostname = gethostname();
cfg.callinfo.pwd      = pwd;

% add information about the function filename and revision to the configuration
cfg.version.name = stack.file;
% the revision number is maintained by SVN in the revision variable in the calling function
if ~exist('revision', 'var')
  cfg.version.id   = 'unknown';
else
  cfg.version.id   = revision;
end

if istrue(ft_getopt(cfg, 'showcallinfo', 'yes'))
  % print some feedback on screen, this is meant to educate the user about
  % the requirements of certain computations and to use that knowledge in
  % distributed computing
  fprintf('the call to "%s" took %d seconds and an estimated %d MB\n', stack.name, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));
end

clear stack
