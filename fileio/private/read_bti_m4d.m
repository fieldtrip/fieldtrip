function [msi] = read_bti_m4d(filename)

% READ_BTI_M4D
%
% Use as
%   msi = read_bti_m4d(filename)

% Copyright (C) 2007, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

[p, f, x] = fileparts(filename);
if ~strcmp(x, '.m4d')
  % add the extension of the header
  filename = [filename '.m4d'];
end

fid = fopen(filename, 'r');
if fid==-1
  error(sprintf('could not open file %s', filename));
end

% start with an empty header structure
msi = struct;

% these header elements contain strings and should be converted in a cell-array
strlist = {
  'MSI.ChannelOrder'
  };

% these header elements contain numbers and should be converted in a numeric array
%   'MSI.ChannelScale'
%   'MSI.ChannelGain'
%   'MSI.FileType'
%   'MSI.TotalChannels'
%   'MSI.TotalEpochs'
%   'MSI.SamplePeriod'
%   'MSI.SampleFrequency'
%   'MSI.FirstLatency'
%   'MSI.SlicesPerEpoch'
% the conversion to numeric arrays is implemented in a general fashion
% and all the fields above are automatically converted
numlist = {};

line = '';

msi.grad.label   = {};
msi.grad.coilpos = zeros(0,3);
msi.grad.coilori = zeros(0,3);

while ischar(line)
  line = cleanline(fgetl(fid));
  if isempty(line) || (length(line)==1 && all(line==-1))
    continue
  end

  sep = strfind(line, ':');
  if length(sep)==1
    key = line(1:(sep-1));
    val = line((sep+1):end);
  elseif length(sep)>1
    % assume that the first separator is the relevant one, and that the
    % next ones are part of the value string (e.g. a channel with a ':' in
    % its name
    sep = sep(1);
    key = line(1:(sep-1));
    val = line((sep+1):end);
  elseif length(sep)<1
    % this is not what I would expect
    error('unexpected content in m4d file');
  end

  if ~isempty(strfind(line, 'Begin')) && (~isempty(strfind(line, 'Meg_Position_Information')) || ~isempty(strfind(line, 'Ref_Position_Information'))) 
    % jansch added the second ~isempty() to accommodate for when the
    % block is about Eeg_Position_Information, which does not pertain to
    % gradiometers, and moreover can be empty (added: Aug 03, 2013)
    
    sep = strfind(key, '.');
    sep = sep(end);
    key = key(1:(sep-1));

    % if the key ends with begin and there is no value, then there is a block
    % of numbers following that relates to the magnetometer/gradiometer information.
    % All lines in that Begin-End block should be treated separately
    val = {};
    lab = {};
    num = {};
    ind = 0;
    while isempty(strfind(line, 'End'))
      line = cleanline(fgetl(fid));
      if isempty(line) || (length(line)==1 && all(line==-1)) || ~isempty(strfind(line, 'End'))
        continue
      end
      ind = ind+1;
      % remember the line itself, and also cut it into pieces
      val{ind} = line;
      % the line is tab-separated and looks like this
      % A68 0.0873437   -0.075789   0.0891512   0.471135    -0.815532   0.336098
      sep = find(line==9); % the ascii value of a tab is 9
      sep = sep(1);
      lab{ind} = line(1:(sep-1));
      num{ind} = str2num(line((sep+1):end));

    end % parsing Begin-End block
    val = val(:);
    lab = lab(:);
    num = num(:);
    num = cell2mat(num);
    
    % the following is FieldTrip specific
    if size(num,2)==6
      msi.grad.label = [msi.grad.label; lab(:)];
      % the numbers represent position and orientation of each magnetometer coil
      msi.grad.coilpos   = [msi.grad.coilpos; num(:,1:3)];
      msi.grad.coilori   = [msi.grad.coilori; num(:,4:6)];
    else
      error('unknown gradiometer design')
    end
  end
  
  % the key looks like 'MSI.fieldname.subfieldname'
  fieldname = key(5:end);

  % remove spaces from the begin and end of the string
  val = strtrim(val);

  % try to convert the value string into something more usefull
  if ~iscell(val)
    % the value can contain a variety of elements, only some of which are decoded here
    if ~isempty(strfind(key, 'Index')) || ~isempty(strfind(key, 'Count')) || any(strcmp(key, numlist))
      % this contains a single number or a comma-separated list of numbers
      val = str2num(val);
    elseif ~isempty(strfind(key, 'Names')) || any(strcmp(key, strlist))
      % this contains a comma-separated list of strings
      val = tokenize(val, ',');
    else
      tmp = str2num(val);
      if ~isempty(tmp)
        val = tmp;
      end
    end
  end

  % assign this header element to the structure
  msi = setsubfield(msi, fieldname, val);

end % while ischar(line)

% each coil weighs with a value of 1 into each channel
msi.grad.tra  = eye(size(msi.grad.coilpos,1));
msi.grad.unit = 'm';

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to remove spaces from the begin and end
% and to remove comments from the lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function line = cleanline(line)
if isempty(line) || (length(line)==1 && all(line==-1))
  return
end
comment = findstr(line, '//');
if ~isempty(comment)
  line(min(comment):end) = ' ';
end
line = strtrim(line);
