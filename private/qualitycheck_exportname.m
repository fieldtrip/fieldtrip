function exportname = exportname_qualitycheck(dataset)

% EXPORTNAME_QUALITYCHECK is used in ft_qualitycheck and the associated cron job

[p, f, x] = fileparts(dataset);

if ft_filetype(dataset, 'ctf_ds')
  try
    info = read_ctf_hist(dataset);
    c = datevec(info.startdate);
    % construct something like 20070320_1230, i.e. YYYYMMDD_HHMM
    exportname = sprintf('%04d%02d%02d_%s%s.mat', c(1), c(2), c(3), info.starttime(1:2), info.starttime(4:5));
  catch
    % the hist file is sometines empty
    % just use the dataset name, but exclude the full path
    exportname = [f, '.mat'];
  end

else
  % just use the dataset name, but exclude the full path
  exportname = [f, '.mat'];
end

