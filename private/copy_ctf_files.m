function copy_ctf_files(oldname, newname, deleteflag)

% COPY_CTF_FILES copies a CTF dataset with all files and directories to a new CTF
% dataset with another name.
%
% Use as
%   copy_brainvision_files(oldname, newname, deleteflag)
%
% Both the old and new name should refer to the CTF dataset directory, including
% the .ds extension.
%
% The third "deleteflag" argument is optional, it should be a boolean
% that specifies whether the original files should be deleted after
% copying or not (default = false).
%
% See also COPY_BRAINVISION_FILES

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

%% deal with inputs

if nargin<3
  deleteflag = false;
end

if ~endsWith(oldname, '.ds') || ~endsWith(newname, '.ds')
  ft_error('you should specify the CTF directory name ending with .ds');
end

if ~isfolder(oldname)
  ft_error('directory "%s" does not exist', oldname);
end

if isfolder(newname)
  ft_error('directory "%s" already exists', oldname);
end

[po, fo, xo] = fileparts(oldname);
[pn, fn, xn] = fileparts(newname);

mkdir(newname);
recurse_copyfile(oldname, newname, fo, fn)
if deleteflag
  rmdir(oldname, 's');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recurse_copyfile(old, new, fo, fn)
list = dir(old);
for i=1:numel(list)
  if strcmp(list(i).name, '.')
    continue
  elseif strcmp(list(i).name, '..')
    continue
  elseif list(i).isdir
    [p, f, x] = fileparts(list(i).name);
    mkdir(fullfile(new, [f, x]));
    if strcmp(f, fo)
      % replace the old directory name by the new one
      recurse_copyfile(fullfile(old, [fo, x]), fullfile(new, [fn, x]), fo, fn);
    else
      % keep the directory name as is
      recurse_copyfile(fullfile(old, [f, x]), fullfile(new, [f, x]), fo, fn);
    end
  else
    [p, f, x] = fileparts(list(i).name);
    if strcmp(f, fo)
      % copy the old file to the new one
      copyfile(fullfile(old, [fo, x]), fullfile(new, [fn, x]));
    elseif ismember(x, {'.cls', '.mrk'})
      % copy the file, search-and-replace the name in the content
      copyfile_and_replace(fullfile(old, [f, x]), fullfile(new, [f, x]), fo, fn)
    else
      % keep the file name and content as is
      copyfile(fullfile(old, [f, x]), fullfile(new, [f, x]));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function copyfile_and_replace(old, new, fo, fn)
fid  = fopen(old, 'r');
f = fread(fid,'*char')';
fclose(fid);
f = strrep(f, fo, fn); % replace the old name by the new one
fid = fopen(new, 'w');
fprintf(fid,'%s',f);
fclose(fid);
