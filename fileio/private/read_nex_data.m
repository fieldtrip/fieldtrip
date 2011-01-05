function [dat] = read_nex_data(filename, hdr, begsample, endsample, chanindx)

% READ_NEX_DATA for Plexon *.nex file
%
% Use as
%   [dat] = read_nex_data(filename, hdr, begsample, endsample, chanindx)
%
% See also READ_NEX_HEADER, READ_NEX_EVENT

% Copyright (C) 2007, Robert Oostenveld
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

try,
  % work with the original header, not the FieldTrip one
  hdr = hdr.orig;
catch
  % assume that we got the original header
end

numsmp = cell2mat({hdr.varheader.numsmp});
adindx = find(cell2mat({hdr.varheader.typ})==5);
smpfrq = hdr.varheader(adindx(1)).wfrequency;
sgn    = chanindx;
nsmp   = (endsample-begsample+1);
dat    = zeros(length(sgn), nsmp);

fid = fopen(filename, 'r', 'ieee-le');
for sgnlop=1:length(sgn)

  if hdr.varheader(sgn(sgnlop)).typ == 0
    % read a spike channel
    status = fseek(fid,hdr.varheader(sgn(sgnlop)).offset,'bof');
    
    % read the sample indices at which spikes occurred
    tim = fread(fid,hdr.varheader(sgn(sgnlop)).cnt,'int32');
    % downsample from 40kHz to the A/D sampling frequency
    tim = (tim ./ hdr.filheader.frequency) * smpfrq + 1;  % add one sample, since ts=0 corresponds to sample=1
    tim = round(tim);   % needed because of the edges in histc
    % select only samples between the desired begin and end
    % tim = tim(find(tim>=begsample & tim<=endsample));
    % convert sample indices into a continuous signal
    if ~isempty(tim)
      dum = histc(tim-begsample, 0:(nsmp-1), 1);
      dat(sgnlop,:) = dum(:)';
    end

  elseif hdr.varheader(sgn(sgnlop)).typ == 5
    % read an A/D channel
    status = fseek(fid,hdr.varheader(sgn(sgnlop)).offset,'bof');
    
    % this just reads the times of LFP starts
    tim = fread(fid,hdr.varheader(sgn(sgnlop)).cnt,'int32');
    % this just reads the indices of LFP starts
    ind = fread(fid,hdr.varheader(sgn(sgnlop)).cnt,'int32');
    if length(ind)>1
      error('multiple A/D segments are not supported');
    end

    % convert from timestamps to samples, expressed in the sampling frequency of the AD channels
    tim = (tim ./ hdr.filheader.frequency) * smpfrq;
    tim = round(tim);
    
    ch_begsample = begsample - tim;
    ch_endsample = endsample - tim;
    
    if (ch_begsample<1)
      error(sprintf('cannot read before the begin of the recorded data (channel %d)', sgn(sgnlop)));
    elseif (ch_endsample>hdr.varheader(sgn(sgnlop)).numsmp)
      error(sprintf('cannot read beyond the end of the recorded data (channel %d)', sgn(sgnlop)));
    end
    
    % seek to the beginning of the interesting data, correct for the A/D card initialisation delay
    fseek(fid,(ch_begsample-1)*2, 'cof');
    % read the actual data for the whole channel
    dum = fread(fid,nsmp,'int16');
    % convert to mV
    dat(sgnlop,:) = dum(:)' * hdr.varheader(sgn(sgnlop)).adtomv;

  else
    % warning('unsupported data format for channel %s', hdr.label{sgn(sgnlop)});
  end

end
status = fclose(fid);

