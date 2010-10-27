function ft_documentationreference(outdir)

% FT_DOCUMENTATIONREFERENCE is a helper function to maintain the online
% documentation

% Copyright (C) 2008, Robert Oostenveld
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
% $Id: ft_preprocessing.m 948 2010-04-21 18:02:21Z roboos $

fieldtripdefs

p = fileparts(which(mfilename));

f1 = dir(fullfile(p, 'ft_*.m'));
f1 = {f1.name}';

f2 = dir(fullfile(p, 'public', '*.m'));
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

funname = cat(1, f1, f2, f3, f4, f5, f6, f7, f8, f9,f10);

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
