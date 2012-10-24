function [dat] = ft_preproc_padding(dat, padtype, padlength)
% FT_PREPROC_PADDING performs padding on the data, i.e. adds or removes 
% samples to/from the data matrix.
%
% Use as
%   [dat] = ft_preproc_padding(dat, padtype, padlength)
% where
%   dat         data matrix (Nchan1 X Ntime)
%   padtype     'zero', 'edge', 'mirror' or 'remove'
%   padlength   scalar, number of samples that will be padded 
%
% 
% See also FT_PREPROCESSING

% Copyright (C) 2012, Jörn M. Horschig, Robert Oostenveld, Jan-Mathijs Schoffelen
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


nchans = size(dat, 1);
nsamples = size(dat, 2);

switch(padtype)  
  case 'remove'
    dat = dat(:, padlength+1:end-padlength);
    return;
  case 'mirror'
    padbeg = 1:padlength;
    padend = 1:padlength;
    begsample = 1;
    endsample = 0;
    while padlength > begsample % this will be a linear piecewise function consisting of two pieces taking turns
      endsample                                 = begsample + min(padlength-endsample, nsamples-1);
      padend(begsample:endsample-1)             = fliplr(mod(0:(endsample-begsample-1), nsamples)+10-endsample+begsample);
      padbeg(end-endsample+2:end-begsample+1)   = fliplr(mod(0:(endsample-begsample-1), nsamples)+2);
      begsample = endsample-1;
      
      if padlength > begsample
        endsample                               = begsample + min(padlength-endsample+1, nsamples-1);
        padend(begsample:endsample)             = mod(0:(endsample-begsample), nsamples)+1;
        padbeg(end-endsample+1:end-begsample+1) = mod(0:(endsample-begsample), nsamples)+10-endsample+begsample;
        begsample = endsample+1;
      end
    end    
    dat = [dat(:, padbeg) dat dat(:, padend)];
    return;
  case 'edge'
    dat = [dat(:,1)*ones(nchans,n) dat dat(:,end)*ones(nchans,n)];
    return;
  case 'zero'
    dat = [zeros(nchans,n) dat zeros(nchans,n)];
    return;
  otherwise
    error('unknown padding option');
end

end