function varargout = cstructdecode(buf, varargin);

% CSTRUCTDECODE decodes a structure from a uint8 buffer
%
% See READ_NEURALYNX_NEV for an example

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: cstructdecode.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.2  2008/02/18 13:36:40  roboos
% added single and double
%
% Revision 1.1  2007/12/18 16:14:38  roboos
% new implementation
%

if ~isa(buf, 'uint8')
  error('incorrect type of input data, should be uint8');
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
      error('incorrect type specification');
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

