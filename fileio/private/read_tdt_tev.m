function [tev, tsq] = read_tdt_tev(filename)

% READ_TDT_TANK
%
% Use as
%   [tev, tsq] = read_tdt_tank(filename)
%
% Note:
%   tev file contains event binary data
%   tev and tsq files work together to get an event's data and attributes
%   sev files contains streamed binary data

% the possible values of the data format are
% DFORM_FLOAT     0
% DFORM_LONG      1
% DFORM_SHORT     2
% DFORM_BYTE      3
% DFORM_DOUBLE    4
% DFORM_QWORD     5

% the possible values of the data type are
% 0         % 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
% 257       % 1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
% 33025     % 1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1
% 33073     % 1     0     0     0     1     1     0     0     1     0     0     0     0     0     0     1
% 34817     % 1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     1
% if bit 16 is set (0x8000) it means that the pointer points into the nev file, otherwise the pointer should be interpreted as value

% the tev file should be accompanied by a tsq file
[p, f, x] = fileparts(filename);
tevfile = filename;
tsqfile = fullfile(p, [f, '.tsq']);

% the tsq file might be accompanied by one or multiple sev files
sevfile = dir(fullfile(p, '*.sev'));
sevfile = {sevfile.name}';

% read the miniblock header information from the tsq file
tsq = read_tdt_tsq(tsqfile);

% this will hold the data snippets
tev = cell(size(tsq));

% open the nev file
fidnev = fopen(tevfile, 'rb');
% open the sev files
fidsev = zeros(size(sevfile));
for i=1:length(sevfile)
  % add the full path
  fidsev(i) = fopen(fullfile(p, sevfile{i}), 'rb');
end

for i=1:numel(tsq)
  if ~bitget(tsq(i).type, 16)
    % the offset does not point into a file, but should be interpreted as value itself
    tev{i} = tsq(i).offset;

  elseif tsq(i).type==33073
    % assume that it points into the sev file

    % this is the part of the sev filename that distinguishes them
    ident = sprintf('%s_ch%d.sev', tsq(i).code, tsq(i).channel);
    % determine whether the corresponding sev file is present
    index = ~cellfun(@isempty, strfind(sevfile, ident));

    if any(index)
      switch tsq(i).format
        case 0
          fmt = 'float32';
        case 1
          fmt = 'int32';
        case 2
          fmt = 'int16';
        case 3
          fmt = 'int8';
        case 4
          fmt = 'double';
        case 5
          ft_error('don''t know what a DFORM_QWORD is');
        otherwise
          ft_error('unknown tsq.type');
      end % switch
      siz = double(tsq(i).size)-10; % in words
      fseek(fidsev(index), tsq(i).offset, 'bof');
      tev{i} = fread(fidsev(index), siz, fmt);
    else
      tev{i} = [];
    end

  else
    % assume that it points into the tev file
    switch tsq(i).format
      case 0
        fmt = 'float32';
      case 1
        fmt = 'int32';
      case 2
        fmt = 'int16';
      case 3
        fmt = 'int8';
      case 4
        fmt = 'double';
      case 5
        ft_error('don''t know what a DFORM_QWORD is');
      otherwise
        ft_error('unknown tsq.type');
    end % switch
    siz = double(tsq(i).size)-10; % in words
    fseek(fidnev, tsq(i).offset, 'bof');
    tev{i} = fread(fidnev, siz, fmt);
  end
end

fclose(fidnev);
for i=1:length(fidsev)
  fclose(fidsev(i));
end

% ucode     = unique({tsq.code});
% uchannel  = unique([tsq.channel]);


