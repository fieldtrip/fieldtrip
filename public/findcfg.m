function [val, status] = findcfg(cfg, var);

% FINDCFG searches for an element in the cfg structure
% or in the nested previous cfgs
%
% Use as
%   [val] = findcfg(cfg, var)
% where the name of the variable should be specified as string.
%
% e.g.
%   trl   = findcfg(cfg, 'trl')
%   event = findcfg(cfg, 'event')

% Copyright (C) 2006, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if var(1)~='.'
  var = ['.' var];
end
val   = [];
depth = 0;
status = 0;

while ~status
  depth = depth + 1;
  if issubfield(cfg,  var)
    val = getsubfield(cfg, var);
    status = 1;
  elseif issubfield(cfg, '.previous');
    [val, status] = findcfg(cfg.previous, var);
     if status, break; end;
  elseif iscell(cfg) 
    for i=1:length(cfg)
      [val, status] = findcfg(cfg{i}, var);
      if status, break; end;
    end
  else
    status = -1;
    break
  end
end

