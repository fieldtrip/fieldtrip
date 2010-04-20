function asc = read_eyelink_asc(filename)

% READ_EYELINK_ASC reads the header information, input triggers, messages
% and all data points from an Eyelink *.asc file
%
% Use as
%   asc = read_eyelink_asc(filename)

% Copyright (C) 2010, Robert Oostenveld

fid = fopen(filename, 'rt');

header  = {};
msg     = {};
input   = [];
dat     = [];
current   = 0;

while ~feof(fid)
  tline = fgetl(fid);

  if regexp(tline, '^[0-9]');
    tmp   = sscanf(tline, '%f');
    nchan = numel(tmp);
    current = current + 1;

    if size(dat,1)<nchan
      % increase the allocated number of channels
      dat(nchan,:) = 0;
    end

    if size(dat, 2)<current
      % increase the allocated number of samples
      dat(:,end+10000) = 0;
    end

    % add the current sample to the data matrix
    dat(1:nchan, current) = tmp;


  elseif regexp(tline, '\*\*.*')
    header = cat(1, header, {tline});


  elseif regexp(tline, '^MSG')
    msg = cat(1, msg, {tline});


  elseif regexp(tline, '^INPUT')
    [val, num] = sscanf(tline, 'INPUT %d %d');
    this.timestamp = val(1);
    this.value     = val(2);
    if isempty(input)
      input = this;
    else
      input = cat(1, input, this);
    end


  else
    % all other lines are not parsed
  end

end

fclose(fid);

asc.header  = header;
asc.msg     = msg;
asc.input   = input;
asc.dat     = dat(:,1:current);
