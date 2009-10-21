function spikefixdmafile(cfg)

% SPIKEFIXDMAFILE fixes the problem in DMA files due to stopping
% and restarting the acquisition. It takes one Neuralynx DMA file and
% and creates seperate DMA files, each corresponding with one continuous
% section of the recording.
%
% Use as
%   spikefixdmafile(cfg);
% where the configuration should contain
%   cfg.dataset   = string with the name of the DMA log file
%   cfg.output    = string with the name of the DMA log file, (default is determined automatic)
%   cfg.numchans  = number of channels (default = 256)
%
% See also SPIKESPLITTING

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: spikefixdmafile.m,v $
% Revision 1.6  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.5  2008/08/13 14:20:59  roboos
% fixed some bugs
%
% Revision 1.4  2008/08/13 10:38:50  roboos
% fixed fprintf bug and number of additional channels to 18
%
% Revision 1.3  2008/08/13 10:22:40  roboos
% changed the logic for jumping in the file
%
% Revision 1.2  2008/08/11 14:11:30  roboos
% added option cfg.numchans
%
% Revision 1.1  2008/06/16 18:47:05  roboos
% new implementation
%

fieldtripdefs

% set the general defaults
if ~isfield(cfg, 'dataset'),  cfg.dataset = [];           end
if ~isfield(cfg, 'output'),   cfg.output = [];            end
if ~isfield(cfg, 'numchans'), cfg.numchans = 256;         end

if isempty(cfg.output)
  [p, f, x] = fileparts(cfg.dataset);
  cfg.output = fullfile(p, f); % without the extension
end

try
  hdr = read_header(cfg.dataset);
catch
  disp(lasterr);
end

nchan     = cfg.numchans+18;
fsample   = 32556;
hdroffset = 16384;
count     = 0;

fin = fopen(cfg.dataset, 'rb');
header = fread(fin, hdroffset, 'uchar=>uchar');


% skip the part that has been processed sofar
fseek(fin, hdroffset, 'bof');

ok = 1;
while (ok)

  % find the beginning of the next section in the file
  buf = zeros(1, nchan);
  while (buf(1)~=2048 && buf(2)~=1)
    % read a single block
    buf = fread(fin, nchan, 'uint32=>uint32');
    if feof(fin) || length(buf)~=nchan
      ok = 0;
      break;
    end
    if (buf(1)~=2048 && buf(2)~=1)
      % jump back to the original location and slide one sample forward to try again
      fseek(fin, -(nchan-1)*4, 'cof');
    else
      % jump back to the original location and continue with the copying to the new file
      fseek(fin, -nchan*4, 'cof');
    end
  end

  % the first file gets an "a" appended, the second a "b", etc.
  count = count + 1;
  filename = sprintf('%s%s.nrd', cfg.output, char(96+count));
  fprintf('writing section %d to "%s"\n', count, filename);
  
  fout = fopen(filename, 'wb');
  fwrite(fout, header, 'uchar');

  while (buf(1)==2048 && buf(2)==1)
    % read a single blockr
    buf = fread(fin, nchan, 'uint32=>uint32');
    if feof(fin) || length(buf)~=nchan
      ok = 0;
      break;
    end
    if (buf(1)==2048 && buf(2)==1)
      % write the block to the other file
      fwrite(fout, buf, 'uint32');
    end
  end

  fclose(fout);

end % while ok

fclose(fin);

