function documentationreference(outdir);

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

f5 = dir(fullfile(p, 'forwinv', '*.m'));
f5 = {f5.name}';

f6 = dir(fullfile(p, 'realtime', 'rt_*.m'));
f6 = {f6.name}';

funname = cat(1, f1, f2, f3, f4, f5, f6);

for i=1:length(funname)
  [p, funname{i}, x] = fileparts(funname{i});
end

% create the desired output directory
if ~isdir(outdir)
mkdir(outdir);
end

for i=1:length(funname)
  filename = fullfile(outdir, [funname{i} '.txt']);
  str = help(funname{i});
  fid = fopen(filename, 'wt');
  fprintf(fid, '=====  %s =====\n\n', upper(funname{i}));
  fprintf(fid, 'Note that this reference documentation is identical to the help that is displayed in Matlab when you type "help %s".\n\n', funname{i});
  fprintf(fid, '<code>\n');  % required for docuwiki
  fprintf(fid, '%s', str);
  fprintf(fid, '</code>\n');  % required for docuwiki
  fclose(fid);
end


