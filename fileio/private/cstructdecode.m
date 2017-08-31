function varargout = cstructdecode(buf, varargin)

% CSTRUCTDECODE decodes a structure from a uint8 buffer
%
% See READ_NEURALYNX_NEV for an example

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

if ~isa(buf, 'uint8')
  ft_error('incorrect type of input data, should be uint8');
end

nbytes = numel(buf);
nfield = length(varargin);

wordsize = zeros(1,nfield);
for i=1:nfield
  switch varargin{i}
  case 'uint8'
    wordsize(i) = 1;
  case 'int8'
    wordsize(i) = 1;
  case 'uint16'
    wordsize(i) = 2;
  case 'int16'
    wordsize(i) = 2;
  case 'uint32'
    wordsize(i) = 4;
  case 'int32'
    wordsize(i) = 4;
  case 'uint64'
    wordsize(i) = 8;
  case 'int64'
    wordsize(i) = 8;
  case {'float32' 'single'}
    varargin{i} = 'single';
    wordsize(i) = 4;
  case {'float64' 'double'}
    varargin{i} = 'double';
    wordsize(i) = 8;
  otherwise
    if strncmp(varargin{i}, 'char', 4)
      if length(varargin{i})>4
        % assume a string like 'char128' which means 128 characters
        wordsize(i) = str2num(varargin{i}(5:end));
        varargin{i} = 'char';
      else
        wordsize(i) = 1;
      end
    else
      ft_error('incorrect type specification');
    end
  end
end

pklen = sum(wordsize);
pknum = nbytes/sum(wordsize);

buf = reshape(buf, pklen, pknum);

for i=1:nfield
  rowbeg = sum(wordsize(1:(i-1)))+1;
  rowend = sum(wordsize(1:(i-0)))+0;
  sel = buf(rowbeg:rowend,:);
  if strcmp(varargin{i}, 'char')
    varargout{i} = char(sel)';
  else
    varargout{i} = typecast(sel(:), varargin{i});
  end
end

