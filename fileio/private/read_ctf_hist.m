function info = read_ctf_hist(filename)

% READ_CTF_HIST

% Copyright (C) 2010, Arjen Stolk

[p, f, x] = fileparts(filename);

if strcmp(x, '.ds')
  % the filename should point to the hist file, not the ds
  filename = fullfile(p, [f x], [f '.hist']);
end

fileline = 0;
fid      = fopen(filename,'r');

while fileline >= 0
  fileline = fgets(fid);
  if ~isempty(findstr(fileline,'Collection started'))
    startdate = sscanf(fileline(findstr(fileline,'Collection started:'):end),'Collection started: %s');
    info.starttime = sscanf(fileline(findstr(fileline,startdate):end),strcat(startdate, '%s'));
    info.startdate = startdate;
  end
  if ~isempty(findstr(fileline,'Collection stopped'))
    stopdate = sscanf(fileline(findstr(fileline,'Collection stopped:'):end),'Collection stopped: %s');
    info.stoptime = sscanf(fileline(findstr(fileline,stopdate):end),strcat(stopdate, '%s'));
    info.stopdate = stopdate;
  end
  if ~isempty(findstr(fileline,'Dataset name'))
    info.datasetname = sscanf(fileline(findstr(fileline,'Dataset name'):end),'Dataset name %s');
  end
end

fclose(fid);
