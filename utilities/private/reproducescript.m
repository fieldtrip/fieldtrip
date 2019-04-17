function reproducescript(reproducescript_dir, filename, highest_ft, cfg, outputcfg)

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

% We're going to explicitly mention any input files going into the pipeline
% at this point in a comment in the script file. We should mention any
% files that (1) have not been produced as an output in the same pipeline
% and that (2) have not been used as an input before. Do this by simply
% checking any .mat file references containing "input" in the stringified cfg against any .mat
% file references already in the script file.
re = ['(?<=' regexptranslate('escape', reproducescript_dir) '[/\\]{1})(\w+_input_\w+(\.mat){1})'];
if exist(filename, 'file')
  script = fileread(filename);
  existing_infiles = regexp(script, re, 'match');
else
  existing_infiles = {};
end
new_infiles = regexp(tmpcfg, re, 'match');
new_infiles = unique(setdiff(new_infiles, existing_infiles));
  
ft_info('writing script to file ''%s''\n', filename);
fid = fopen(filename, 'a+');
fprintf(fid, '%%%%\n\n');

for k = 1:numel(new_infiles)
  fprintf(fid, '%% a new input variable is entering the pipeline here: %s\n', new_infiles{k});
end
if numel(new_infiles) > 0
  fprintf(fid, '\n');
end
  
fprintf(fid, 'cfg = [];\n');
fprintf(fid, '%s\n', tmpcfg);

if outputcfg
  % this is for ft_definetrial, ft_artifact_zvalue, etc.
  fprintf(fid, 'cfg = %s(cfg);\n\n', highest_ft);
else
  fprintf(fid, '%s(cfg);\n\n', highest_ft);
end

fclose(fid);
