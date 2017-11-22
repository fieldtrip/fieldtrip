function [dat] = ft_preproc_padding(dat, padtype, prepadlength, postpadlength)

% FT_PREPROC_PADDING performs padding on the data, i.e. adds or removes samples
% to or from the data matrix.
%
% Use as
%   [dat] = ft_preproc_padding(dat, padtype, padlength)
% or as
%   [dat] = ft_preproc_padding(dat, padtype, prepadlength, postpadlength)
% where
%   dat           data matrix (Nchan x Ntime)
%   padtype       'zero', 'mean', 'localmean', 'edge', 'mirror', 'nan' or 'remove'
%   padlength     scalar, number of samples that will be padded
%   prepadlength  scalar, number of samples that will be padded before the data
%   postpadlength scalar, number of samples that will be padded after the data
%
% If padlength is used instead of prepadlength and postpadlength, padding
% will be symmetrical (i.e. padlength = prepadlength = postpadlength)
%
% See also FT_PREPROCESSING

% Copyright (C) 2012, J?rn M. Horschig, Robert Oostenveld, Jan-Mathijs Schoffelen
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

if nargin < 4
  postpadlength = prepadlength;
end

if prepadlength == 0 && postpadlength == 0
  return;
end

nchans = size(dat, 1);
nsamples = size(dat, 2);

switch(padtype)
  case 'remove'
    dat = dat(:, prepadlength+1:end-postpadlength);
    
  case 'mirror'
    padbeg = 1:prepadlength;
    padend = 1:postpadlength;
    
    % predata padding
    begsample = 1;
    endsample = 0;
    while prepadlength > begsample % this will be a linear piecewise function consisting of two pieces taking turns
      endsample                                 = begsample + min(prepadlength-endsample, nsamples-1);
      padbeg(end-endsample+2:end-begsample+1)   = fliplr(mod(0:(endsample-begsample-1), nsamples)+2);
      begsample = endsample-1;
      
      if prepadlength > begsample
        endsample                               = begsample + min(prepadlength-endsample+1, nsamples-1);
        padbeg(end-endsample+1:end-begsample+1) = mod(0:(endsample-begsample), nsamples)+nsamples-endsample+begsample;
        begsample = endsample+1;
      end
    end
    
    % postdata padding
    begsample = 1;
    endsample = 0;
    while postpadlength > begsample % this will be a linear piecewise function consisting of two pieces taking turns
      endsample                                 = begsample + min(postpadlength-endsample, nsamples-1);
      padend(begsample:endsample-1)             = fliplr(mod(0:(endsample-begsample-1), nsamples)+nsamples-endsample+begsample);
      begsample = endsample-1;
      
      if postpadlength > begsample
        endsample                               = begsample + min(postpadlength-endsample+1, nsamples-1);
        padend(begsample:endsample)             = mod(0:(endsample-begsample), nsamples)+1;
        begsample = endsample+1;
      end
    end
    
    dat       = [dat(:, padbeg) dat dat(:, padend)];
    
  case 'edge'
    dat       = [dat(:,1)*ones(1,prepadlength) dat dat(:,end)*ones(1,postpadlength)];
    
  case 'mean'
    mu        = mean(dat, 2);
    dat       = [mu*ones(1,prepadlength) dat mu*ones(1,postpadlength)];
    
  case 'localmean'
    prepad    = min(prepadlength, floor(size(dat, 2)/2));
    edgeleft  = mean(dat(:, 1:prepad), 2);
    postpad   = min(postpadlength, floor(size(dat, 2)/2));
    edgeright = mean(dat(:, 1+end-postpad:end), 2);
    dat       = [edgeleft*ones(1,prepadlength) dat edgeright*ones(1,postpadlength)];
    
  case 'zero'
    dat       = [zeros(nchans,prepadlength) dat zeros(nchans,postpadlength)];
    
  case 'nan'
    dat       = [nan(nchans,prepadlength) dat nan(nchans,postpadlength)];
    
  otherwise
    ft_error('unknown padding option');
end

end

