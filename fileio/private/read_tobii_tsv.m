function out = read_tobii_tsv(filename, varargin)

% READ_TOBII_TSV
%
% Use as
%   hdr = read_tobii_tsv(filename)
% or
%   dat = read_tobii_tsv(filename, tsv, begsample, endsample)

needhdr = (nargin==1);
needdat = (nargin>1);


if needhdr
  % get the original header details
  
  fid = fopen(filename, 'rt');
  tsv = struct();
  ln  = 0;
  while ln<30 && ~feof(fid)
    ln   = ln+1;
    line = fgetl(fid);
    if any(line==':') && line(end)~=':'
      % this looks like a header line
      [key, val] = strtok(line, ':');
      key = fixname(key);
      val = val(2:end); % skip the ':'
      val = deblank2(val);
      tsv.(key) = val;
    end
    
    if ~isempty(strfind(line, 'Timestamp')) && ~isempty(strfind(line, 'Gaze'))
      % this looks like the line with the headings of all the columns
      tsv.header = tokenize(line, '\t');
      % remember where the data starts
      tsv.headerline  = ln;
      tsv.datapointer = ftell(fid);
    end
    
  end
  
  fclose(fid);
  
  hdr = [];
  hdr.orig = tsv;
  
  % return the FieldTrip header information
  out = hdr;
  
elseif needdat
  
  hdr       = varargin{1};
  begsample = varargin{2};
  endsample = varargin{3};
  
  % get the original header details
  tsv = hdr.orig;
  fid = fopen(filename, 'rt');
  fseek(fid, tsv.datapointer, 'bof');
  
  cursample = 1;
  while cursample<begsample
    line = fgetl(fid);
    curline = curline + 1;
  end
  dat = zeros(length(tsv.header), endsample-begsample+1);
  while cursample<endsample
    line = fgetl(fid);
    val = tokenize(line, 9); % horizontal tab is 9 in ascii table
    num = cellfun(@str2num, val, 'UniformOutput', false);
    valid = cellfun(@isempty, num, 'UniformOutput', true);
    val(~str) = num;
    val( str) = nan;
    dat(:,cursample-begsample+1) = val;
    curline = curline + 1;
  end
  fclose(fid);
  
  % return the data
  out = hdr;
end

