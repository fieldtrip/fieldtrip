function savevar(filename, varname, value)

% LOADVAR is a helper function for cfg.inputfile

% Copyright (C) 2010, Robert Oostenveld
%
% $Id$

fprintf('writing ''%s'' to file ''%s''\n', varname, filename);

eval(sprintf('%s = value;', varname));
save(filename, varname);
