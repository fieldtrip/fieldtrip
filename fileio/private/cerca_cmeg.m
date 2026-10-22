function varargout = cerca_cmeg(filename, hdr, begsample, endsample, chanindx)

% CERCA_CMEG Read a Cerca OPM-MEG binary file and associated metadata.
%
% Use as
%   hdr = cerca_cmeg(filename);
%   dat = cerca_cmeg(filename, hdr, begsample, endsample, chanindx);
%   evt = cerca_cmeg(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT
% See also BIDS_TSV, BIOPAC_ACQ, BUCN_TXT, EEGSYNTH_TSV, EVENTS_TSV, LIBERTY_CSV, MAUS_TEXTGRID, MOTION_C3D, OPENBCI_TXT, OPENPOSE_KEYPOINTS, OPENSIGNALS_TXT, OPENVIBE_MAT, OPM_FIL, QUALISYS_TSV, QUSPIN_LVM, SCCN_XDF, SENSYS_CSV, SNIRF, SPIKEGLX_BIN, UNICORN_CSV, XSENS_MVNX

% Copyright (C) 2026 Robert Oostenveld
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

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

% use the full filename including path to distinguish between similarly named files in different directories
filename = which(filename);

Adim_conv = [2^32, 2^16, 2^8, 1];

nbytes = dir(filename);
nbytes = nbytes.bytes;

% ---- Read first 8 bytes as block header -------------------------
fid = fopen(filename, 'r', 'ieee-be');
if fid < 0
  error('Cannot open %s', filename);
end
blockhdr = fread(fid, 8, 'uint8=>double');
blockhdr = reshape(blockhdr, [4, 2])';    % 2x4, big-endian order
blockdim = Adim_conv * blockhdr';         % 1x2 : [rows, cols]
nelements = 1 + prod(blockdim);
bytesperblock = 8 * nelements;
nblocks = nbytes / bytesperblock;

% ---- Magic number ------------------------------------------------
frewind(fid);
magic = fread(fid, 1, 'double=>double', 0, 'ieee-be');

if needhdr
  % ---- Locate auxiliary files --------------------------------------
  [p, f, x] = fileparts(filename);
  basename  = fullfile(p, f(1:end-4)); % remove _meg

  megjson         = [basename '_meg.json'];
  channelstsv     = [basename '_channels.tsv'];
  sessioninfo     = [basename '_SessionInfo.txt'];
  helmconfig      = [basename '_HelmConfig.tsv'];
  sensortransform = [basename '_SensorTransform.tsv'];

  if isfile(megjson)
    megjson = ft_read_json(megjson);
  else
    fprintf('cannot find %s\n', megjson);
    megjson = [];
  end

  if isfile(channelstsv)
    ws = warning;
    warning off MATLAB:table:ModifiedAndSavedVarnames
    channelstsv = ft_read_tsv(channelstsv);
    for i=1:numel(channelstsv.name)
      % remove spaces and brackets
      channelstsv.name{i}(channelstsv.name{i}==' ') = [];
      channelstsv.name{i}(channelstsv.name{i}=='[') = [];
      channelstsv.name{i}(channelstsv.name{i}==']') = [];
      % remove the trailing Z from non-OPM channels
      if startsWith(channelstsv.name{i}, 'Trigger') && endsWith(channelstsv.name{i}, 'Z')
        channelstsv.name{i}(end) = [];
      elseif startsWith(channelstsv.name{i}, 'BNC') && endsWith(channelstsv.name{i}, 'Z')
        channelstsv.name{i}(end) = [];
      end
    end
    warning(ws)
  else
    fprintf('cannot find %s\n', channelstsv);
    channelstsv = [];
  end

  if isfile(sessioninfo)
    fprintf('reading ''%s''\n', sessioninfo);
    fid = fopen(sessioninfo, 'rt');
    sessioninfo = fread(fid, [1 inf], 'char=>char');
    fclose(fid);
  else
    fprintf('cannot find %s\n', sessioninfo);
    sessioninfo = [];
  end

  if isfile(helmconfig)
    tsv = ft_read_tsv(helmconfig);
    tsv = tsv(1:end-1,:); % drop the last line

    helmetlab  = tsv{:,1};   % labels for the helmet slots
    pos        = tsv{:,2:4};
    ori        = tsv{:,5:7};
    xy         = tsv{:,8:9};
    channellab = tsv{:,10};  % labels for the OPM sensor channels

    % remove spaces from the helmet slot and OPM channel names
    for i=1:numel(helmetlab)
      helmetlab{i}(isspace(helmetlab{i})) = [];
    end
    for i=1:numel(channellab)
      channellab{i}(isspace(channellab{i})) = [];
    end

    % select only the relevant rows from the helmet configuration TSV file
    chansel = ~strcmp(channellab, 'none');
    helmetlab = helmetlab(chansel);
    channellab = channellab(chansel);
    pos = pos(chansel,:);
    ori = ori(chansel,:);
    xy  = xy(chansel,:);

    % construct the grad so that it can be in the output
    grad = [];
    grad.label    = channellab;
    grad.coilpos  = pos;
    grad.coilori  = ori;
    grad.tra      = eye(length(grad.label));
    grad.coordsys = 'RAS';
    grad.unit     = 'm';

    % also store the whole table
    helmconfig = tsv;
  else
    fprintf('cannot find %s\n', helmconfig);
    grad = [];
    helmconfig = [];
  end

  if isfile(sensortransform)
    fprintf('reading ''%s''\n', sensortransform);
    fid = fopen(sensortransform, 'rt');
    sensortransform = fscanf(fid, '%f', [4 4]);
    fclose(fid);
  else
    fprintf('cannot find %s\n', sensortransform);
    sensortransform = [];
  end

  % construct the FieldTrip header
  hdr = [];
  hdr.label       = [{'Timestamp'}; channelstsv.name];
  hdr.chantype    = [{'unknown'}; channelstsv.type];
  hdr.chanunit    = [{'unknown'}; channelstsv.units];
  hdr.Fs          = megjson.SamplingFrequency;
  hdr.nChans      = blockdim(1);
  hdr.nSamples    = blockdim(2) * nblocks;
  hdr.nTrials     = 1; % continuous
  hdr.nSamplesPre = 0;

  if ismember('V_nT', channelstsv.Properties.VariableNames)
    converted = ~isnan(channelstsv.V_nT);
    converted = [false; converted]; % the timestamp channel is not converted
    for i=find(converted')
      % the data from this channel will be converted from V to T in the code below
      hdr.chanunit{i} = 'T';
    end
  else
    warning('cannot convert raw OPM voltages (V) to magnetic field strength (T)')
  end

  if ~isempty(grad)
    hdr.grad = grad;
  end

  % also store the original metadata
  hdr.orig.megjson         = megjson;
  hdr.orig.channelstsv     = channelstsv;
  hdr.orig.sessioninfo     = sessioninfo;
  hdr.orig.helmconfig      = helmconfig;
  hdr.orig.sensortransform = sensortransform;

  varargout = {hdr};

elseif needdat

  % ---- Read data blocks --------------------------------------------
  dat = zeros(blockdim(1), blockdim(2), floor(nblocks));
  for ii = 1:floor(nblocks)
    fseek(fid, (ii-1)*bytesperblock, 'bof');
    tmp = fread(fid, nelements, 'double=>double', 0, 'ieee-be');
    assert(tmp(1) == magic, 'Magic number mismatch in block %d', ii);
    tmp = tmp(2:end);
    tmp = reshape(tmp, fliplr(blockdim));
    dat(:, :, ii) = tmp';
  end
  fclose(fid);

  % reshape to [nrows, ncols*nblocks] in Fortran order
  dat = reshape(dat, blockdim(1), []);

  if ismember('V_nT', hdr.orig.channelstsv.Properties.VariableNames)
    % convert from V to T, this is also reflected in the hdr.chanunit above
    scale = hdr.orig.channelstsv.V_nT;
    scale(isnan(scale)) = scale(isnan(scale)) * 1e9; % convert from nT to T
    scale = [1; scale]; % for the timestamp channel
    scale(isnan(scale)) = 1; % for the non-OPM channels
    scale = sparse(diag(1./scale));
    dat = scale * dat;
  end

  % make the selection of channels and samples
  dat = dat(chanindx, begsample:endsample);

  varargout = {dat};

elseif needevt
  % detect the flanks in the trigger channels
  evt = read_trigger(filename, 'header', hdr, 'chanindx', find(strcmpi(hdr.chantype, 'TRIG')));

  varargout = {evt};
end
