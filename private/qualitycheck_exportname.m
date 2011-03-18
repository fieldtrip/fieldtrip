function exportname = exportname_qualitycheck(dataset)

% EXPORTNAME_QUALITYCHECK is used in ft_qualitycheck and the associated cron job

[p, f, x] = fileparts(dataset);
histfile = fullfile(p, [f x], [f '.hist']); % this only applies to some CTF datasets

if isdir(dataset) && exist(histfile, 'file')
  info = read_ctf_hist(histfile);
  c = datevec(info.startdate);
  % construct something like 20070320_1230, i.e. YYYYMMDD_HHMM
  exportname = sprintf('%04d%02d%02d_%s%s.mat', c(1), c(2), c(3), info.starttime(1:2), info.starttime(4:5));
else
  % just use the dataset name, but exclude the full path
  exportname = [f, '.mat'];
end

