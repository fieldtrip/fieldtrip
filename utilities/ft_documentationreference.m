function ft_documentationreference(outdir)

% FT_DOCUMENTATIONREFERENCE is a function to maintain the online reference
% documentation.
%
% Normal users will not be calling this function, but will rather look at
% http://www.fieldtriptoolbox.org/reference where the output of this function can
% be found.
%
% See also FT_DOCUMENTATIONCONFIGURATION, MATLAB2MARKDOWN

% Copyright (C) 2008-2019, Robert Oostenveld
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

ft_defaults

[ftver, ftpath] = ft_version;

subdir = {
  '.'
  'connectivity'
  'contrib/misc'
  'contrib/nutmegtrip'
  'contrib/spike'
  'engine'
  'external/artinis'
  'fileio'
  'forward'
  'inverse'
  'peer'
  'plotting'
  'preproc'
  'qsub'
  'realtime/example'
  'realtime/online_eeg'
  'realtime/online_meg'
  'realtime/online_mri'
  'specest'
  'statfun'
  'trialfun'
  'utilities'
  };

infile = {};

% find all functions that should be included in the reference documentation
for i=1:length(subdir)
  f = dir(fullfile(ftpath, subdir{i}, '*.m'));
  f = {f.name}';
  for j=1:length(f)
    f{j} = fullfile(ftpath, subdir{i}, f{j});
  end
  infile = cat(1, infile, f);
end

% create the desired output directory
if ~isfolder(outdir)
  mkdir(outdir);
end

for i=1:length(infile)
  [p, f, x] = fileparts(infile{i});
  outfile = fullfile(outdir, [f '.md']);
  matlab2markdown(infile{i}, outfile, 'monospacehelp', true, 'overwrite', true, 'pageheader', 'jekyll', 'pagetitle', f, 'highlight', 'plaintext');
end
