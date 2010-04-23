function update_pcode

funlist = {
  'addCTFtrial'
  'getCTFBalanceCoefs'
  'getCTFdata'
  'readCPersist'
  'readCTFMRI'
  'readCTFds'
  'readCTFhdm'
  'setCTFDataBalance'
  'writeCPersist'
  'writeCTFMRI'
  'writeCTFds'
  'writeCTFhdm'
  'writeMarkerFile'
  'writeRes4'
  };

outdir = 'pcode';

% create the desired output directory
if ~isdir(outdir)
  mkdir(outdir);
end

% create the m-files with opnly the help in the output dir
for i=1:length(funlist)
  filename = fullfile(outdir, [funlist{i} '.m']);
  str = help(funlist{i});
  fid = fopen(filename, 'wt');
  % FIXME does not fully work
  fprintf(fid, '%s', str);
  fclose(fid);
end

% create the p-files
dir = pwd;
cd(outdir);
for i=1:length(funlist)
  filename = fullfile(dir, [funlist{i} '.m']);
  pcode(filename);
end
cd(dir);