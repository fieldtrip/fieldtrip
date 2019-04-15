function reproducescript(filename, cfg, outputcfg)

% This is a helper function to create a script that reproduces the analysis. It
% appends the configuration and the function call to a MATLAB script.

% Copyright (C) 2019, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if isfield(cfg, 'callinfo') && isfield(cfg.callinfo, 'usercfg')
  tmpcfg = removefields(cfg.callinfo.usercfg, ignorefields('reproducescript'));
else
  tmpcfg = removefields(cfg, ignorefields('reproducescript'));
end

tmpcfg = copyfields(cfg, tmpcfg, {'inputfile', 'outputfile'});
tmpcfg = printstruct('cfg', tmpcfg);

ft_info('writing script to file ''%s''\n', filename);

fid = fopen(filename, 'a+');

fprintf(fid, '%%%%\n\n');
fprintf(fid, 'cfg = [];\n');
fprintf(fid, '%s\n', tmpcfg);

st = dbstack(2);

if outputcfg
  % this is for ft_definetrial, ft_artifact_zvalue, etc.
  fprintf(fid, 'cfg = %s(cfg);\n\n', st(2).name);
else
  fprintf(fid, '%s(cfg);\n\n', st(2).name);
end

fclose(fid);
