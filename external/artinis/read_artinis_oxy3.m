function data = read_artinis_oxy3(filename, header, begsample, endsample, chanindx)
% reads Artinix oxy3-files into fieldtrip format
% use as
%   hdr   = read_artinis_oxy3(filename)
% or 
%   event = read_artinis_oxy3(filename, read_event)
% where read_event is a boolean. If 'true', the function returns events. 
% If 'false' the function returns the header. 
% or 
%   data  = read_artinis_oxy3(filename, header, [begsample], [endsample], [chanindx])
% where begsample, endsample and chanindx are optional. The returned
% variables will be in FieldTrip style. 
%
% See also FT_READ_HEADER, FT_READ_DATA

% Copyright (c) 2015 by Artinis Medical Systems BV
% Author    Jörn M. Horschig
% E-Mail    askforinfo@artinis.com
% Web       www.artinis.com

if nargin > 2 && islogical(header)
  % header must be defined if begsample and endsample are used
  error('wrong type of header variable.')
end

if nargin == 1 || nargin == 2 && islogical(header) && ~header    
  data = read_oxy3_header(filename);
elseif nargin == 2  && islogical(header) && header
  data = read_oxy3_event(filename);
else % nargin > 1 && ~islogical(header)
  data = read_oxy3_data(filename, header);
  
  if nargin <5
    chanindx = 1:size(data, 1);
    if nargin < 4
      endsample = size(data, 2);    
      if nargin < 3
        begsample = 1;
      end
    end
  end
        
  data = data(chanindx, begsample:endsample);
end