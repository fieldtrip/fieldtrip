function E = read_buffer_offline_events(eventfile, hdr)
% function E = read_buffer_offline_events(eventfile, header)
%
% This function reads FCDC buffer-type events from a binary file.

% (C) 2010 S. Klanke

type = {
  'char'
  'uint8'
  'uint16'
  'uint32'
  'uint64'
  'int8'
  'int16'
  'int32'
  'int64'
  'single'
  'double'
};

wordsize = {
  1 % 'char'
  1 % 'uint8'
  2 % 'uint16'
  4 % 'uint32'
  8 % 'uint64'
  1 % 'int8'
  2 % 'int16'
  4 % 'int32'
  8 % 'int64'
  4 % 'single'
  8 % 'double'
};

F = fopen(eventfile,'rb',hdr.orig.endianness);
numEvt = 0;

fseek(F, 0, 'eof');
fileSize = ftell(F);
fseek(F, 0, 'bof');

maxEvt = floor(fileSize / 32);
if maxEvt == 0
   E = [];
   return;
end

E = repmat(struct('type', 0, 'sample', 0, 'duration', 0, 'offset', 0, 'value', 0), [maxEvt,1]);

num=0;

while ~feof(F)
  try
    type_type   = fread(F, 1, 'uint32');
    type_numel  = fread(F, 1, 'uint32');
    value_type  = fread(F, 1, 'uint32');
    value_numel = fread(F, 1, 'uint32');
    sample      = fread(F, 1, 'int32');
    offset      = fread(F, 1, 'int32');
    duration    = fread(F, 1, 'int32');
    bufsize     = fread(F, 1, 'uint32');
    buf         = fread(F, bufsize, 'uint8=>uint8');
  catch
    break;
  end
  
  type_size  = wordsize{type_type+1}*type_numel;
  value_size = wordsize{value_type+1}*value_numel;
  
  if type_size + value_size > bufsize
     warning('Invalid event definition - skipping remaining contents');
     break;
  end
    
  num=num+1;
  E(num).sample = sample + 1; % offset are 1-based in Matlab
  E(num).offset = offset;
  E(num).duration = duration;
  if type_type==0
    E(num).type = char(buf(1:type_size))';
  else
    E(num).type = typecast(buf(1:type_size), type{type_type+1});
  end
  if value_type==0
    E(num).value = char(buf((type_size+1):(type_size+value_size)))';
  else
    E(num).value = typecast(buf((type_size+1):(type_size+value_size)), type{value_type+1});
  end
end

E = E(1:num);
