% FT_POSTAMBLE_CALLINFO

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
cfg.callinfo.pwd      = pwd;

% add information about the function filename and revision to the configuration
cfg.version.name = stack.file;
cfg.version.id   = revision;

% give some feedback
fprintf('the call to "%s" took %d seconds and an estimated %d MB\n', stack.name, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));

clear stack