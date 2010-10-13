function savevar(filename, varname, value)

% SAVEVAR is a helper function for cfg.outputfile

% Copyright (C) 2010, Robert Oostenveld
%
% $Id$

fprintf('writing ''%s'' to file ''%s''\n', varname, filename);

eval(sprintf('%s = value;', varname));
save(filename, varname, '-v7.3');
