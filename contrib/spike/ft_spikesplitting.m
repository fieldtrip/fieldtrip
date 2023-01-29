function [cfg] = ft_spikesplitting(cfg)

% FT_SPIKESPLITTING reads a single Neuralynx DMA log file and writes each
% individual channel to a separate file.
%
% Use as
%   [cfg] = ft_spikesplitting(cfg)
%
% The configuration should contain
%   cfg.dataset   = string with the name of the DMA log file
%   cfg.channel   = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.output    = string with the name of the splitted DMA dataset directory, (default is determined automatic)
%   cfg.latency   = [begin end], (default = 'all')
%   cfg.feedback  = string, (default = 'textbar')
%   cfg.format    = 'int16' or 'int32' (default = 'int32')
%   cfg.downscale = single number or vector (for each channel), corresponding to the number of bits removed from the LSB side (default = 0)
%
% This function expects the DMA file to be read as AD units (and not in uV)
% and will write the same AD values to the splitted DMA files. If you
% subsequently want to process the splitted DMA, you should look up the
% details of the headstage amplification and the Neuralynx amplifier and
% scale the values accordingly.
%
% See also FT_SPIKEDOWNSAMPLE, FT_SPIKEANALYSIS

% Copyright (C) 2007-2008, Robert Oostenveld
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


% set the general defaults
if ~isfield(cfg, 'dataset'),          cfg.dataset = [];                 end
if ~isfield(cfg, 'channel'),          cfg.channel = 'all';              end
if ~isfield(cfg, 'latency'),          cfg.latency = [0 inf];            end
if ~isfield(cfg, 'feedback'),         cfg.feedback = 'textbar';         end
if ~isfield(cfg, 'output'),           cfg.output = [];                  end % see below
if ~isfield(cfg, 'format'),           cfg.format = 'int32';             end
if ~isfield(cfg, 'downscale'),        cfg.downscale = 0;                end
if ~isfield(cfg, 'headerformat'),     cfg.headerformat = [];            end

if isempty(cfg.output)
  % set smart defaults for the output
  [p, f, x]  = fileparts(cfg.dataset);
  cfg.output = fullfile(p, [f '.sdma']);
end

status = mkdir(cfg.output);
if ~status
  error('error creating splitted DMA output dataset %s', cfg.output);
end
fprintf('writing to output directory ''%s''\n', cfg.output);

% read the header of the completete dataset
hdr = ft_read_header(cfg.dataset, 'headerformat', cfg.headerformat);

if isfield(hdr, 'orig') && isfield(hdr.orig, 'Header')
  [p, f, x]  = fileparts(cfg.output);
  headerfile = fullfile(cfg.output, [f '.txt']);
  fprintf('writing header information to %s\n', headerfile);
  fid = fopen(headerfile, 'wb');
  fwrite(fid, hdr.orig.Header, 'char');
  fclose(fid);
end

if hdr.nSamples<1
  error('the input dataset contains no samples');
end

if numel(cfg.downscale)==1
  % each channel should have its own downscale factor
  cfg.downscale = repmat(cfg.downscale, 1, hdr.nChans);
end

% determine the selected channels
cfg.channel        = ft_channelselection(cfg.channel, hdr.label);
chansel            = match_str(hdr.label, cfg.channel);  % this is a list with the indices of the selected channels
writechan          = zeros(1,hdr.nChans);
writechan(chansel) = 1;                                  % this is a logical/boolean vector with a 0/1 for each channel

if numel(cfg.downscale) ~= numel(writechan)
  error('the downscale factor should be specified for all channels, not only for the selected ones');
end

if isequal(cfg.latency, 'all')
  begsample       = 1;
  endsample       = hdr.nSamples;
elseif numel(cfg.latency)==2
  begsample       = max(round(cfg.latency(1) * hdr.Fs + 1), 1);
  endsample       = min(round(cfg.latency(2) * hdr.Fs    ), hdr.nSamples);
else
  errror('incorrect specification of cfg.latency');
end

numsample       = endsample - begsample + 1;
cfg.latency(1)  = (begsample-1)/hdr.Fs;
cfg.latency(2)  = (endsample  )/hdr.Fs;

% give some feedback, based on the complete data
fprintf('data contains %10d samples\n', hdr.nSamples);
fprintf('data selected %10d samples\n', numsample);

s = floor(hdr.nSamples ./ hdr.Fs);
m = floor(s/60);
h = floor(m/60);
m = m - 60*h;
s = s - 60*m - 60*60*h;
fprintf('duration of data      %02dh:%02dm:%02ds\n', h, m, s);

s = floor(numsample ./ hdr.Fs);
m = floor(s/60);
h = floor(m/60);
m = m - 60*h;
s = s - 60*m - 60*60*h;
fprintf('duration of selection %02dh:%02dm:%02ds\n', h, m, s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% split the DMA file into separate channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~ft_filetype(cfg.dataset, 'neuralynx_dma')
  error('unsupported data format for DMA splitting');
end

statuslabel = {
  'stx'
  'pid'
  'siz'
  'tsh'
  'tsl'
  'cpu'
  'ttl'
  'x01'
  'x02'
  'x03'
  'x04'
  'x05'
  'x06'
  'x07'
  'x08'
  'x09'
  'x10'
  'crc'
  };

% the output files will have an extension that is based on the channel name like "dataset.channame.bin"
[p, f, x] = fileparts(cfg.output);
filename  = cell(hdr.nChans,1);
format    = cell(hdr.nChans,1);
for i=1:hdr.nChans
  filename{i} = fullfile(cfg.output, [f '.' hdr.label{i} '.bin']);
  if strmatch(hdr.label{i}, statuslabel)
    format{i} = 'uint32';
  elseif strcmp(cfg.format, 'int16')
    format{i} = 'int16';
  elseif strcmp(cfg.format, 'int32')
    format{i} = 'int32';
  else
    error('incorrect specification of cfg.format');
  end
end

% the status channels are written in their original format
for i=1:hdr.nChans
  if strcmp(format{i}, 'uint32') && cfg.downscale(i)~=0
    warning('not downscaling status channel (%s)', hdr.label{i});
    cfg.downscale(i) = 0;
  end
end

% read and write the data in one-second segments
begsample = max(cfg.latency(1) * hdr.Fs + 1, 1);
endsample = min(cfg.latency(2) * hdr.Fs + 1, hdr.nSamples);
segment   = begsample:hdr.Fs:endsample;
if segment(end)~=endsample
  segment(end+1) = endsample;  % the last segment should end on the end sample
end

% The first version of this file format contained in the first 8 bytes the
% channel label as string. Subsequently it contained 32 bit integer values.
%
% The second version of this file format starts with 8 bytes describing (as
% a space-padded string) the data type. The channel label is contained in
% the filename as dataset.chanlabel.bin.
%
% The third version of this file format starts with 7 bytes describing (as
% a zero-padded string) the data type, followed by the 8th byte which
% describes the downscaling for the 8 and 16 bit integer representations.
% The downscaling itself is represented as uint8 and should be interpreted as
% the number of bits to shift. The channel label is contained in the
% filename as dataset.chanlabel.bin.

% open the output files, one for each selected channel
fid = zeros(hdr.nChans,1);
for j=1:hdr.nChans
  if writechan(j)
    fid(j) = fopen(filename{j}, 'wb', 'ieee-le');   % open the file
    magic = format{j};                              % this used to be the channel name
    magic((end+1):8) = 0;                           % pad with zeros
    magic(8) = cfg.downscale(j);                    % number of bits to shift
    fwrite(fid(j), magic(1:8));                     % write the 8-byte file header
  end
end

ft_progress('init', cfg.feedback, 'splitting data');
for i=1:(length(segment)-1)
  % read one segment of data
  begsample = segment(i);
  endsample = segment(i+1)-1;  % the begin of the next segment minus one
  ft_progress(i/(length(segment)-1), 'splitting data segment %d from %d\n', i, length(segment)-1);
  buf = read_neuralynx_dma(cfg.dataset, begsample, endsample, 'all');
  if ~isa(buf, 'int32')
    error('the buffer is expected to be int32');
  end
  % apply the scaling, this corresponds to bit shifting
  for j=1:hdr.nChans
    buf(j,:) = buf(j,:) ./ (2^cfg.downscale(j));
  end
  % write the segment of data to each of the 274 output files
  % the data is casted into the desired output format
  for j=1:hdr.nChans
    if writechan(j)
      if strcmp(format{j}, 'uint32')
        % the individual file header specifies that it should be interpreted as uint32
        % write the int32 data from the buffer as it is, during reading the sign bit will again be correctly interpreted
        fwrite(fid(j), buf(j,:), 'int32', 'ieee-le');
      else
        fwrite(fid(j), buf(j,:), format{j}, 'ieee-le');
      end
    end
  end
end
ft_progress('close');

% close all output files
for j=1:hdr.nChans
  if writechan(j)
    fclose(fid(j));
  end
end

% do the general cleanup and bookkeeping at the end of the function

ft_postamble provenance

