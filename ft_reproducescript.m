function ft_reproducescript(cfg)

% FT_REPRODUCESCRIPT is a helper function to clean up the script and intermediate
% datafiles that are the result from using the cfg.reproducescript option. You should
% call this function all the way at the end of your analysis. This function will look
% at all intermediate files in the output directory, remove input and output files
% that are the same and update the script accordingly.
%
% Use as
%   ft_reproducescript(cfg)
%
% The configuration structure should contain
%   cfg.reproducescript = string, directory with the script and intermediate data
%
% See also FT_ANALYSISPIPELINE, FT_DEFAULTS

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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init

% get the filename of the script
f = dir(fullfile(cfg.rewritescript, 'script.m'));
assert(numel(f)==1)
f = fullfile(f.folder, f.name);
fid = fopen(f, 'r');
str = fread(fid, inf, 'char=>char')';
fclose(fid);

% get a list of data files
d = dir(fullfile(cfg.rewritescript, '*.mat'));
% get the short names and the full names
ds = {d.name};
df = cellfun(@fullfile, {d.folder}, {d.name}, 'UniformOutput', false);

md5 = cellfun(@CalcMD5, df, 'UniformOutput', false);
for i=1:numel(md5)
  md5{i} = 'a';
end

[c, ia, ic] = unique(md5, 'stable');

% update all file names in the script
for i=1:numel(d)
  str = strrep(str, ds{i}, ds{ic(i)});
end

fid = fopen(f, 'w');
fprintf(fid, '%s', str);
fclose(fid);

