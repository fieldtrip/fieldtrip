function documentationreference(outdir);

fieldtripdefs

p = fileparts(which(mfilename));

f1 = dir(fullfile(p, 'ft_*.m'));
f1 = {f1.name}';

funname = cat(1, f1);

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


