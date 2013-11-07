function ft_documentationreference(outdir)

% FT_DOCUMENTATIONREFERENCE is a function to maintain the online documentation.
%
% Normal users will not be calling this function, but will rather look at
% http://fieldtrip.fcdonders.nl/reference where the output of this function can
% be found.
%
% See FT_DOCUMENTATIONINDEX

% Copyright (C) 2008-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

p = fileparts(which(mfilename));

f1 = dir(fullfile(p, '*.m'));
f1 = {f1.name}';

f2 = dir(fullfile(p, 'utilities', '*.m'));
f2 = {f2.name}';

f3 = dir(fullfile(p, 'preproc', '*.m'));
f3 = {f3.name}';

f4 = dir(fullfile(p, 'fileio', '*.m'));
f4 = {f4.name}';

f5 = dir(fullfile(p, 'forward', '*.m'));
f5 = {f5.name}';

f6 = dir(fullfile(p, 'inverse', '*.m'));
f6 = {f6.name}';

f7 = dir(fullfile(p, 'realtime', '*.m'));
f7 = {f7.name}';

f8 = dir(fullfile(p, 'realtime', 'datasource', '*.m'));
f8 = {f8.name}';

f9 = dir(fullfile(p, 'peer', '*.m'));
f9 = {f9.name}';

f10 = dir(fullfile(p, 'plotting', '*.m'));
f10 = {f10.name}';

f11 = dir(fullfile(p, 'statfun', '*.m'));
f11 = {f11.name}';

f12 = dir(fullfile(p, 'specest', '*.m'));
f12 = {f12.name}';

f13 = dir(fullfile(p, 'connectivity', '*.m'));
f13 = {f13.name}';

f14 = dir(fullfile(p, 'contrib', '*.m'));
f14 = {f14.name}';

f15 = dir(fullfile(p, 'qsub', '*.m'));
f15 = {f15.name}';

f16 = dir(fullfile(p, 'contrib', 'spike', '*.m'));
f16 = {f16.name}';

f17 = dir(fullfile(p, 'contrib', 'trentool', '*.m'));
f17 = {f17.name}';

f18 = dir(fullfile(p, 'utility', '*.m'));
f18 = {f18.name}';

f19 = dir(fullfile(p, 'engine', '*.m'));
f19 = {f19.name}';

funname = cat(1, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f19);

for i=1:length(funname)
  [p, funname{i}, x] = fileparts(funname{i});
end

% create the desired output directory
if ~isdir(outdir)
  mkdir(outdir);
end

funname = flipdim(funname,1); % to avoid problems with overlapping function names
for i=1:length(funname)
  filename = fullfile(outdir, [funname{i} '.txt']);
  str = help(funname{i});
  
  % make text html-compatible
  str = strrep(str, '<', '&lt;');
  str = strrep(str, '>', '&gt;');
  
  % add crossrefs
  for f=1:length(funname)
    str = strrep(str, [' ', upper(funname{f})], [' <a href=/reference/', funname{f}, '><font color=green>', upper(funname{f}),'</font></a>']);
  end
  
  fid = fopen(filename, 'wt');
  fprintf(fid, '=====  %s =====\n\n', upper(funname{i}));
  fprintf(fid, 'Note that this reference documentation is identical to the help that is displayed in Matlab when you type "help %s".\n\n', funname{i});
  fprintf(fid, '<html><pre>\n');   % required for docuwiki > use html preformatted style
  fprintf(fid, '%s', str);
  fprintf(fid, '</pre></html>\n'); % required for docuwiki
  fclose(fid);
end
