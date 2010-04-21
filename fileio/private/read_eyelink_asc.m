function asc = read_eyelink_asc(filename)

% READ_EYELINK_ASC reads the header information, input triggers, messages
% and all data points from an Eyelink *.asc file
%
% Use as
%   asc = read_eyelink_asc(filename)

% Copyright (C) 2010, Robert Oostenveld

fid = fopen(filename, 'rt');

asc.header  = {};
asc.msg     = {};
asc.input   = [];
asc.sfix    = {};
asc.efix    = {};
asc.ssacc   = {};
asc.esacc   = {};
asc.dat     = [];
current   = 0;

while ~feof(fid)
  tline = fgetl(fid);

  if regexp(tline, '^[0-9]');
    tmp   = sscanf(tline, '%f');
    nchan = numel(tmp);
    current = current + 1;

    if size(asc.dat,1)<nchan
      % increase the allocated number of channels
      asc.dat(nchan,:) = 0;
    end

    if size(asc.dat, 2)<current
      % increase the allocated number of samples
      asc.dat(:,end+10000) = 0;
    end

    % add the current sample to the data matrix
    asc.dat(1:nchan, current) = tmp;


  elseif regexp(tline, '^INPUT')
    [val, num] = sscanf(tline, 'INPUT %d %d');
    this.timestamp = val(1);
    this.value     = val(2);
    if isempty(asc.input)
      asc.input = this;
    else
      asc.input = cat(1, asc.input, this);
    end


  elseif regexp(tline, '\*\*.*')
    asc.header = cat(1, asc.header, {tline});


  elseif regexp(tline, '^MSG')
    asc.msg = cat(1, asc.msg, {tline});


  elseif regexp(tline, '^SFIX')
    asc.sfix = cat(1, asc.sfix, {tline});


  elseif regexp(tline, '^EFIX')
    asc.efix = cat(1, asc.efix, {tline});


  elseif regexp(tline, '^SSACC')
    asc.ssacc = cat(1, asc.ssacc, {tline});


  elseif regexp(tline, '^ESACC')
    asc.esacc = cat(1, asc.esacc, {tline});


  else
    % all other lines are not parsed
  end

end

% close the file?
fclose(fid);

% remove the samples that were not filled with real data
asc.dat = asc.dat(:,1:current);
