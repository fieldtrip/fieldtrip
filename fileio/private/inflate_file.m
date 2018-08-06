function outputfile = inflate_file(inputfile)

% INFLATE_FILE helper function to uncompress a compressed file of arbitrary
% compression type. Returns the full path to the extracted file or
% directory, which will be located in a temporary location.

% Copyright (C) 2012, Eelke Spaak
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

% use a cache to prevent extracting the same file multiple times
persistent extractcache

inputfilemd5 = fixname(CalcMD5(inputfile));
if isfield(extractcache, inputfilemd5)
  outputfile = extractcache.(inputfilemd5);
  if exist(outputfile, 'file') % checks for files or folders
    fprintf('compressed dataset has already been extracted to %s\n', outputfile);
    return;
  else
    % apparently the file is gone
    extractcache = rmfield(extractcache, inputfilemd5);
  end
end

% determine compression type
if filetype_check_extension(inputfile, 'zip')
  type = 'zip';
elseif filetype_check_extension(inputfile, 'tar')     ...
    || filetype_check_extension(inputfile, '.tar.gz') ...
    || filetype_check_extension(inputfile, 'tgz')
  type = 'tar';
elseif filetype_check_extension(inputfile, 'gz')
  type = 'gzip';
else
  ft_error('unsupported compression type, only zip/gz/tar/tgz/tar.gz are supported');
end

% determine temporary output folder
outputdir = [tempdir() inputfilemd5];

% give some feedback
fprintf('extracting compressed dataset to %s...\n', [outputdir filesep]);

% do appropriate inflation
switch (type)
  case 'zip'
    unzip(inputfile, outputdir);
  case 'tar'
    untar(inputfile, outputdir);
  case 'gzip'
    gunzip(inputfile, outputdir);
end

% if we extracted only a single file (or directory), return path to that file, if it was a
% set of files, return path to the biggest file
outfiles = dir(outputdir);
if numel(outfiles) == 3 % first two will always be . and ..
  outputfile = [outputdir filesep outfiles(3).name];
else
  siz = [outfiles.bytes];
  maxind = find(siz == max(siz),1); % use find() because we want at most 1
  outputfile = [outputdir filesep outfiles(maxind).name];
end

fprintf('extracted dataset is located at %s\n', outputfile);
extractcache.(inputfilemd5) = outputfile;
