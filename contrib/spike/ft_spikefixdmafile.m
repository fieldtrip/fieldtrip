function ft_spikefixdmafile(cfg)

% FT_SPIKEFIXDMAFILE fixes the problem in DMA files due to stopping and
% restarting the acquisition. It takes one Neuralynx DMA file and and
% creates seperate DMA files, each corresponding with one continuous
% section of the recording.
%
% Use as
%   ft_spikefixdmafile(cfg)
% where the configuration should contain
%   cfg.dataset   = string with the name of the DMA log file
%   cfg.output    = string with the name of the DMA log file, (default is determined automatic)
%   cfg.numchans  = number of channels (default = 256)
%
% See also FT_SPIKESPLITTING

% Copyright (C) 2008, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig

% set the general defaults
if ~isfield(cfg, 'dataset'),  cfg.dataset = [];           end
if ~isfield(cfg, 'output'),   cfg.output = [];            end
if ~isfield(cfg, 'numchans'), cfg.numchans = 256;         end
if ~isfield(cfg, 'headerformat'), cfg.headerformat = [];  end

if isempty(cfg.output)
  [p, f, x] = fileparts(cfg.dataset);
  cfg.output = fullfile(p, f); % without the extension
end

try
  hdr = ft_read_header(cfg.dataset, 'headerformat', cfg.headerformat);
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

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble provenance

