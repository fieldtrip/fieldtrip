function inspect_realtime_jitter_data

filename = 'buffer://localhost:1972';

%% write data to the buffer, you should have the buffer running

% the code below emulates a regular heartbeat on channel one
% with a beat slightly faster than every second

hdr = [];
hdr.Fs            = 300;
hdr.nSamples      = 0;
hdr.nSamplesPre   = 0;
hdr.label         = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '11', '12', '13', '14', '15', '16'};
hdr.nChans        = length(hdr.label);

blocksize = 0.1;                     % in seconds
blocksize = round(blocksize*hdr.Fs); % in samples

block = 0;

while true
  block     = block + 1;
  begsample = block*blocksize + 1;
  endsample = begsample + blocksize - 1;
  dat       = rand(hdr.nChans, blocksize);
  dat(1,:)  = mod(10+(begsample:endsample), 280)==0;
  
  ft_write_data(filename, dat, 'header', hdr, 'append', block~=1);
  
  disp(block);
  pause(blocksize/hdr.Fs); % in seconds
end
