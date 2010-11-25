function tsq = read_tdt_tsq(filename, begblock, endblock)

% READ_TDT_TSQ reads the headers from a Tucker_Davis_technologies TSQ file
%
% tsq file is a heap of event headers, which is ?40 byte each,
% ordered strictly by time
%
% Use as
%   tsq = read_tdt_tsq(filename, begblock, endblock)

% Copyright (C) 2010, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% Here are the possible values of the data format long:
% DFORM_FLOAT     0
% DFORM_LONG      1
% DFORM_SHORT     2
% DFORM_BYTE      3
% DFORM_DOUBLE    4
% DFORM_QWORD     5

if nargin<2 || isempty(begblock)
  begblock = 1;
end

if nargin<3 || isempty(endblock)
  endblock = inf;
end

fid = fopen(filename, 'rb');
offset = (begblock-1)*40; % bytes
fseek(fid, offset, 'cof');
buf = fread(fid, [40, (endblock-begblock+1)], 'uint8=>uint8');
fclose(fid);

buf_size      = buf(1:4,:);
buf_type      = buf(5:8,:);
buf_code      = buf(9:12,:);
buf_channel   = buf(13:14,:);
buf_sortcode  = buf(15:16,:);
buf_timestamp = buf(17:24,:);
buf_offset    = buf(25:32,:);
buf_format    = buf(33:36,:);
buf_frequency = buf(37:40,:);

tsq_size       = typecast(buf_size     (:), 'uint32');
tsq_type       = typecast(buf_type     (:), 'uint32');
tsq_code       = char(buf_code');  % this can sometimes be interpreted as 4 chars
tsq_channel    = typecast(buf_channel  (:), 'uint16');
tsq_sortcode   = typecast(buf_sortcode (:), 'uint16');
tsq_timestamp  = typecast(buf_timestamp(:), 'double');
tsq_offset     = typecast(buf_offset   (:), 'uint64');  % or double
tsq_format     = typecast(buf_format   (:), 'uint32');
tsq_frequency  = typecast(buf_frequency(:), 'single');

tsq = struct(  'size'      , num2cell(tsq_size      ), ... 
               'type'      , num2cell(tsq_type      ), ... 
               'code'      , cellstr(tsq_code      ), ... 
               'channel'   , num2cell(tsq_channel   ), ... 
               'sortcode'  , num2cell(tsq_sortcode  ), ... 
               'timestamp' , num2cell(tsq_timestamp ), ... 
               'offset'    , num2cell(tsq_offset    ), ... 
               'format'    , num2cell(tsq_format    ), ... 
               'frequency' , num2cell(tsq_frequency ));

