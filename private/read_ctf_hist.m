function hist = read_ctf_hist(filename)

% READ_CTF_HIST

% Copyright (C) 2010, Arjen Stolk

[p, f, x] = fileparts(filename);

if strcmp(x, '.ds')
  % the filename should point to the hist file, not the ds
  filename = fullfile(p, [f x], [f '.hist']);
end

fileline = 0;
fid      = fopen(filename,'r');

% the hist file is empty sometimes, in which case the actual reading does not return any useful information
hist = [];
hist.starttime   = 'unknown';
hist.startdate   = 'unknown';
hist.stoptime    = 'unknown';
hist.stopdate    = 'unknown';
hist.datasetname = 'unknown';

while fileline >= 0
  fileline = fgets(fid);
  if ~isempty(findstr(fileline,'Collection started'))
    startdate = sscanf(fileline(findstr(fileline,'Collection started:'):end),'Collection started: %s');
    hist.starttime = sscanf(fileline(findstr(fileline,startdate):end),strcat(startdate, '%s'));
    hist.startdate = startdate;
  end
  if ~isempty(findstr(fileline,'Collection stopped'))
    stopdate = sscanf(fileline(findstr(fileline,'Collection stopped:'):end),'Collection stopped: %s');
    hist.stoptime = sscanf(fileline(findstr(fileline,stopdate):end),strcat(stopdate, '%s'));
    hist.stopdate = stopdate;
  end
  if ~isempty(findstr(fileline,'Dataset name'))
    hist.datasetname = sscanf(fileline(findstr(fileline,'Dataset name'):end),'Dataset name %s');
  end
end

fclose(fid);
